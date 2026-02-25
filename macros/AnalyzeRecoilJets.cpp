// AnalyzeRecoilJets.cpp
//
// Usage (ROOT):
//   root -l -q AnalyzeRecoilJets.cpp
// -----------------------------------------------------------------------------

#include "AnalyzeRecoilJets.h"

namespace ARJ
{
  // =============================================================================
  // Analysis modules (implementation lives here; shared utilities live in the header)
  // =============================================================================
  namespace analysis
  {
    void RunEventLevelQA(Dataset& ds)
      {
        cout << ANSI_BOLD_CYN << "\n==============================\n"
             << "[SECTION 1] EVENT-LEVEL + INVENTORY (" << ds.label << ")\n"
             << "==============================" << ANSI_RESET << "\n";

        const double nEvt = ReadEventCount(ds);
        cout << ANSI_BOLD_RED
             << "[" << ds.label << "] accepted events = " << std::fixed << std::setprecision(0) << nEvt
             << ANSI_RESET << "\n";

        vector<InvItem> items;
        CollectInventoryRecursive(ds.topDir, ds.topDirName + "/", items);
        std::sort(items.begin(), items.end(), [](const InvItem& a, const InvItem& b){ return a.path < b.path; });
        PrintInventoryToTerminal(ds, items);

        TH1* hV = GetObj<TH1>(ds, "h_vertexZ", true, true, true);
        if (!hV)
        {
          cout << ANSI_BOLD_YEL << "[WARN] Missing h_vertexZ for " << ds.label << ANSI_RESET << "\n";
          return;
        }

        TH1F* hFixed = RebinToFixedBinWidthVertexZ(hV, vzCutCm);
        if (!hFixed)
        {
          cout << ANSI_BOLD_YEL << "[WARN] Failed to build fixed-binwidth vertexZ for " << ds.label << ANSI_RESET << "\n";
          return;
        }

        string outPath;
        if (ds.isSim)
          outPath = JoinPath(ds.outBase, "GeneralEventLevelQA/zvtx_SIM.png");
        else
          outPath = JoinPath(ds.outBase, "baselineData/GeneralEventLevelQA/zvtx_DATA_" + ds.trigger + ".png");

        vector<string> lines;
        lines.push_back(TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data());
        DrawAndSaveTH1_Common(ds, hFixed, outPath, "v_{z} [cm]", "Counts", lines, false);

        // ------------------------------------------------------------------
        // SIM ONLY: truth-vs-reco vertex 2D QA (pre-cut fill in RecoilJets.cc)
        //   X = v_z^truth, Y = v_z^(reco-used)
        // ------------------------------------------------------------------
        if (ds.isSim)
        {
          TH2* h2 = GetObj<TH2>(ds, "h_vzTruthVsReco", true, true, true);
          if (!h2)
          {
            cout << ANSI_BOLD_YEL << "[WARN] Missing h_vzTruthVsReco for " << ds.label << ANSI_RESET << "\n";
          }
          else
          {
            const string out2 = JoinPath(ds.outBase, "GeneralEventLevelQA/vzTruthVsReco_SIM.png");

            vector<string> l2;
            l2.push_back("Pre-cut fill (before |v_{z}| cut)");
            l2.push_back("X: v_{z}^{truth}  |  Y: v_{z}^{reco-used}");

            DrawAndSaveTH2_Common(ds, h2, out2,
                                 "v_{z}^{truth} [cm]",
                                 "v_{z}^{reco-used} [cm]",
                                 "Counts",
                                 l2,
                                 false);
          }
        }

        delete hFixed;
      }

      // =============================================================================
      // Trigger turn-on QA (pp DATA only): doNotScale hMaxClusterEnergy overlays + ratios
      // =============================================================================
      void RunTriggerAna_DoNotScaleMaxClusterEnergy(Dataset& ds)
      {
        if (ds.isSim) return;

        static std::set<std::string> s_doneOutDirs;

        const std::string ppBase  = DirnameFromPath(ds.outBase);
        const std::string outDir  = JoinPath(ppBase, "triggerQA");

        if (s_doneOutDirs.count(outDir)) return;
        s_doneOutDirs.insert(outDir);

        EnsureDir(outDir);

        const std::string mbdShort = "MBD_NandS_geq_1";
        const std::string p3Short  = "Photon_3_GeV_plus_MBD_NS_geq_1";
        const std::string p4Short  = "Photon_4_GeV_plus_MBD_NS_geq_1";

        const std::string prefix  = "h_maxEnergyClus_NewTriggerFilling_doNotScale_";

        auto getHist = [&](const std::string& trigShort)->TH1*
        {
          if (!ds.file) return nullptr;

          TDirectory* dir = ds.file->GetDirectory(trigShort.c_str());
          if (!dir) return nullptr;

          const std::string hname = prefix + trigShort;
          TH1* h = dynamic_cast<TH1*>(dir->Get(hname.c_str()));
          return h;
        };

        TH1* hMBD = getHist(mbdShort);
        TH1* hP3  = getHist(p3Short);
        TH1* hP4  = getHist(p4Short);

        if (!hMBD || !hP3 || !hP4)
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] Missing doNotScale trigger hist(s) needed for triggerAna outputs.\n"
               << "       Need directories + hists:\n"
               << "         " << mbdShort << "/" << prefix << mbdShort << "\n"
               << "         " << p3Short  << "/" << prefix << p3Short  << "\n"
               << "         " << p4Short  << "/" << prefix << p4Short  << "\n"
               << ANSI_RESET << "\n";
          return;
        }

        // ------------------------------------------------------------------
        // (1) Overlay: MBD vs Photon3 vs Photon4 (doNotScale max cluster energy)
        // ------------------------------------------------------------------
        {
          TCanvas c("c_trigAna_overlay", "c_trigAna_overlay", 900, 700);
          c.cd();

          hMBD->SetTitle("");
          hMBD->GetXaxis()->SetTitle("Cluster Energy [GeV]");
          hMBD->GetYaxis()->SetTitle("Entries");

          hMBD->SetLineWidth(2);
          hP3->SetLineWidth(2);
          hP4->SetLineWidth(2);

          hMBD->Draw("HIST");
          hP3->Draw("HIST SAME");
          hP4->Draw("HIST SAME");

          TLegend leg(0.55, 0.68, 0.88, 0.88);
          leg.SetBorderSize(0);
          leg.SetFillStyle(0);
          leg.AddEntry(hMBD, "MBD N&S >= 1", "l");
          leg.AddEntry(hP3,  "Photon 3 GeV + MBD NS >= 1", "l");
          leg.AddEntry(hP4,  "Photon 4 GeV + MBD NS >= 1", "l");
          leg.Draw();

          const std::string outPng = JoinPath(outDir, "hMaxClusterEnergy_doNotScale_overlay.png");
          c.SaveAs(outPng.c_str());
          cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";
        }

        // ------------------------------------------------------------------
        // (2) Ratios: Photon3/MBD and Photon4/MBD (turn-on style)
        // ------------------------------------------------------------------
        {
          TH1* rP3 = dynamic_cast<TH1*>(hP3->Clone("ratio_P3_over_MBD"));
          TH1* rP4 = dynamic_cast<TH1*>(hP4->Clone("ratio_P4_over_MBD"));

          if (!rP3 || !rP4)
          {
            cout << ANSI_BOLD_YEL << "[WARN] Failed to clone histograms for ratio plot." << ANSI_RESET << "\n";
            return;
          }

          rP3->SetDirectory(nullptr);
          rP4->SetDirectory(nullptr);

          rP3->Divide(hP3, hMBD, 1.0, 1.0, "B");
          rP4->Divide(hP4, hMBD, 1.0, 1.0, "B");

          TCanvas c("c_trigAna_ratio", "c_trigAna_ratio", 900, 700);
          c.cd();

          rP3->SetTitle("");
          rP3->GetXaxis()->SetTitle("Cluster Energy [GeV]");
          rP3->GetYaxis()->SetTitle("Efficiency");

          rP3->SetLineWidth(2);
          rP4->SetLineWidth(2);

          rP3->SetMinimum(0.0);
          rP3->SetMaximum(1.2);

          rP3->Draw("HIST");
          rP4->Draw("HIST SAME");

          TLegend leg(0.55, 0.72, 0.88, 0.88);
          leg.SetBorderSize(0);
          leg.SetFillStyle(0);
          leg.AddEntry(rP3, "Photon 3 GeV + MBD NS >= 1  /  MBD", "l");
          leg.AddEntry(rP4, "Photon 4 GeV + MBD NS >= 1  /  MBD", "l");
          leg.Draw();

          const std::string outPng = JoinPath(outDir, "hMaxClusterEnergy_doNotScale_ratioToMBD.png");
          c.SaveAs(outPng.c_str());
          cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";

          delete rP3;
          delete rP4;
        }
      }

      // =============================================================================
      // Section 2: preselection fail counters (terminal only)
      // =============================================================================
      void RunPreselectionFailureTable(Dataset& ds)
      {
          cout << ANSI_BOLD_CYN << "\n==============================\n"
               << "[SECTION 2] PRESELECTION FAIL QA (TERMINAL + PLOTS) (" << ds.label << ")\n"
               << "==============================" << ANSI_RESET << "\n";

          // ---------------------------------------------------------------------------
          // Output directory for preselection plots:
          //   SIM  -> <outBase>/PurityABCD/preselection/
          //   DATA -> <outBase>/baselineData/PurityABCD/preselection/
          //          (where <outBase> is already the per-trigger folder)
          // ---------------------------------------------------------------------------
          string outDir;
          if (ds.isSim) outDir = JoinPath(ds.outBase, "PurityABCD/preselection");
          else          outDir = JoinPath(ds.outBase, "baselineData/PurityABCD/preselection");
          EnsureDir(outDir);

          const int wBin = 10;
          const int wN   = 12;

          cout << std::left
               << std::setw(wBin) << "pTbin"
               << std::right
               << std::setw(wN) << "wetaFail"
               << std::setw(wN) << "et1Low"
               << std::setw(wN) << "et1High"
               << std::setw(wN) << "et1Out"
               << std::setw(wN) << "e11e33Hi"
               << std::setw(wN) << "e32e35Lo"
               << std::setw(wN) << "e32e35Hi"
               << std::setw(wN) << "e32e35Out"
               << "\n";
          cout << string(wBin + 8*wN, '-') << "\n";

          // Requested x-axis labels (keep these exact)
          const char* xLabels[8] = {
            "#frac{E_{11}}{E_{33}}",
            "et1 < 0.6",
            "et1 > 1.0",
            "0.6 < et1 < 1.0",
            "#frac{E_{32}}{E_{35}} < 0.8",
            "#frac{E_{32}}{E_{35}} > 1.0",
            "0.8 < #frac{E_{32}}{E_{35}} < 1.0",
            "w_{#eta}^{cogX} < 0.6"
          };

          // Clean, distinct solid colors per bin
          const int binColors[8] = {
            kGray+1,
            kGreen+2,
            kRed+1,
            kMagenta+1,
            kOrange+7,
            kYellow+2,
            kAzure+1,
            kCyan+1
          };

          // Keep values so we can also build a 3x3 summary table at the end.
          vector< vector<double> > valsByPt(kNPtBins, vector<double>(8, 0.0));

          // ---------------------------------------------------------------------------
          // Helper: draw one preselection-bar plot into the *current pad* (gPad).
          // If keepAlive != nullptr, created histograms are pushed there and MUST be
          // deleted by the caller AFTER the final SaveCanvas.
          // ---------------------------------------------------------------------------
          auto DrawPreselectionBarsIntoPad =
            [&](const PtBin& b, const double vals[8], bool compact, vector<TObject*>* keepAlive)
          {
            if (!gPad) return;

            // pad layout
            if (compact)
            {
                gPad->SetLeftMargin(0.17);
                gPad->SetRightMargin(0.08);
                gPad->SetTopMargin(0.14);
                gPad->SetBottomMargin(0.46);
            }
            else
            {
              gPad->SetLeftMargin(0.12);
              gPad->SetRightMargin(0.05);
              gPad->SetTopMargin(0.12);
              gPad->SetBottomMargin(0.34);
            }
            gPad->SetTicks(1,1);

            double ymax = 0.0;
            for (int ib = 0; ib < 8; ++ib) ymax = std::max(ymax, vals[ib]);
            const double yMaxPlot = (ymax > 0.0) ? (1.35 * ymax) : 1.0;

            // Axis-only histogram (frame + labels)
            TH1F* hAxis = new TH1F(
              TString::Format("h_preFailAxis_%s_%s_%s",
                ds.label.c_str(), b.folder.c_str(), compact ? "tbl" : "one").Data(),
              "",
              8, 0.5, 8.5
            );
            hAxis->SetDirectory(nullptr);
            hAxis->SetStats(0);
            hAxis->SetMinimum(0.0);
            hAxis->SetMaximum(yMaxPlot);

            for (int ib = 1; ib <= 8; ++ib) hAxis->GetXaxis()->SetBinLabel(ib, xLabels[ib-1]);

            hAxis->GetYaxis()->SetTitle("Fail counts");
            hAxis->GetXaxis()->SetTitle("");
            hAxis->GetXaxis()->LabelsOption("v");

            hAxis->GetXaxis()->SetLabelSize(compact ? 0.065 : 0.055);
            hAxis->GetXaxis()->SetLabelOffset(compact ? 0.018 : 0.012);

            // Y-axis title: bigger + closer in the 3x3 table
            hAxis->GetYaxis()->SetTitleSize(compact ? 0.065 : 0.055);
            hAxis->GetYaxis()->SetTitleOffset(compact ? 1.05 : 1.05);

            hAxis->SetLineColor(1);
            hAxis->SetLineWidth(2);
            hAxis->SetFillStyle(0);
            hAxis->Draw("hist");

            // Draw solid 2D bars (NO BAR2 -> avoids 3D shading)
            // Keep alive until canvas is saved (pad stores pointers).
            for (int ib = 1; ib <= 8; ++ib)
            {
              TH1F* hb = new TH1F(
                TString::Format("h_preFailBar_%s_%s_%s_b%d",
                  ds.label.c_str(), b.folder.c_str(), compact ? "tbl" : "one", ib).Data(),
                "",
                8, 0.5, 8.5
              );
              hb->SetDirectory(nullptr);
              hb->SetStats(0);

              hb->SetBinContent(ib, vals[ib-1]);

              hb->SetFillStyle(1001);
              hb->SetFillColor(binColors[ib-1]);
              hb->SetLineColor(1);
              hb->SetLineWidth(2);

              hb->SetBarWidth(0.90);
              hb->SetBarOffset(0.05);

              hb->Draw("BAR SAME");

              if (keepAlive) keepAlive->push_back(hb);
              else delete hb; // only safe if caller saves immediately after draw
            }

            // Numeric counts above each bar
            TLatex t;
            t.SetTextFont(42);
            t.SetTextAlign(22); // centered
            t.SetTextSize(compact ? 0.055 : 0.034);


            for (int ib = 1; ib <= 8; ++ib)
            {
              const double y = vals[ib-1];
              if (y <= 0.0) continue;

              const double x = hAxis->GetXaxis()->GetBinCenter(ib);
              const double yText = std::min(y + 0.03*yMaxPlot, 0.95*yMaxPlot);
              t.DrawLatex(x, yText, TString::Format("%.0f", y).Data());
            }

            // Clean top-left annotation
            vector<string> box;
            if (!compact)
            {
              vector<string> hdr = DefaultHeaderLines(ds);
              if (!hdr.empty()) box.push_back(hdr[0]);
              box.push_back("Preselection fails (inclusive)");
              box.push_back(TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
              box.push_back("NOTE: one photon can increment multiple categories");
              DrawLatexLines(0.15, 0.93, box, 0.032, 0.040);
            }
            else
            {
              // --- Compact table annotation (NDC) ---
              TLatex tt;
              tt.SetTextFont(42);
              tt.SetNDC();

              // pT range (top-left) — larger, replaces where "inclusive fails" used to sit
              tt.SetTextAlign(13);   // left, top
              tt.SetTextSize(0.07);
              tt.DrawLatex(0.2, 0.8,
                TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data()
              );

              // "Inclusive Fails" (top-right)
              tt.SetTextAlign(33);   // right, top
              tt.SetTextSize(0.062);
              tt.DrawLatex(0.88, 0.93, "Inclusive Fails (Preselection)");
            }

            if (keepAlive) keepAlive->push_back(hAxis);
            else delete hAxis;
          };

          const auto& bins = PtBins();

          // ---------------------------------------------------------------------------
          // Per pT bin: print terminal row + write an individual PNG
          // ---------------------------------------------------------------------------
          for (int i = 0; i < kNPtBins; ++i)
          {
            const PtBin& b   = bins[i];
            const string suf = b.suffix;

            const double weta = Read1BinCount(ds, "h_preFail_weta" + suf);
            const double et1L = Read1BinCount(ds, "h_preFail_et1_low" + suf);
            const double et1H = Read1BinCount(ds, "h_preFail_et1_high" + suf);
            const double et1O = et1L + et1H;

            const double e11H = Read1BinCount(ds, "h_preFail_e11e33_high" + suf);
            const double e32L = Read1BinCount(ds, "h_preFail_e32e35_low" + suf);
            const double e32H = Read1BinCount(ds, "h_preFail_e32e35_high" + suf);
            const double e32O = e32L + e32H;

            // ---- terminal row (UNCHANGED) ----
            cout << std::left << std::setw(wBin) << b.label
                 << std::right
                 << std::setw(wN) << std::fixed << std::setprecision(0) << weta
                 << std::setw(wN) << et1L
                 << std::setw(wN) << et1H
                 << std::setw(wN) << et1O
                 << std::setw(wN) << e11H
                 << std::setw(wN) << e32L
                 << std::setw(wN) << e32H
                 << std::setw(wN) << e32O
                 << "\n";

            // Values in the exact requested order
            valsByPt[i][0] = e11H;
            valsByPt[i][1] = et1L;
            valsByPt[i][2] = et1H;
            valsByPt[i][3] = et1O;
            valsByPt[i][4] = e32L;
            valsByPt[i][5] = e32H;
            valsByPt[i][6] = e32O;
            valsByPt[i][7] = weta;

            // Individual plot canvas
            TCanvas c(
              TString::Format("c_preFail_%s_%s", ds.label.c_str(), b.folder.c_str()).Data(),
              "c_preFail", 1300, 720
            );

            c.cd();

            // IMPORTANT: keep objects alive until after SaveCanvas
            vector<TObject*> keep;
            {
              double vv[8];
              for (int ib = 0; ib < 8; ++ib) vv[ib] = valsByPt[i][ib];
              DrawPreselectionBarsIntoPad(b, vv, false, &keep);
            }

            const string outPng = JoinPath(
              outDir,
              TString::Format("preselectionFails_bar_%s.png", b.folder.c_str()).Data()
            );
            SaveCanvas(c, outPng);

            // cleanup objects used in this canvas
            for (auto* obj : keep) delete obj;
          }

          // ---------------------------------------------------------------------------
          // 3x3 summary table (pT bins in order)
          // ---------------------------------------------------------------------------
          {
            TCanvas cTbl(
                TString::Format("c_preFail_table3x3_%s", ds.label.c_str()).Data(),
                "c_preFail_table3x3", 2400, 1400
            );
            cTbl.Divide(3,2, 0.002, 0.002);

            vector<TObject*> keepTbl;
            for (int i = 0; i < kNPtBins; ++i)
            {
              cTbl.cd(i+1);

              double vv[8];
              for (int ib = 0; ib < 8; ++ib) vv[ib] = valsByPt[i][ib];
              DrawPreselectionBarsIntoPad(bins[i], vv, true, &keepTbl);
            }

            // Save one 3x3 PNG
            SaveCanvas(cTbl, JoinPath(outDir, "table3x3_preselectionFails.png"));

            // cleanup objects used in the table canvas
            for (auto* obj : keepTbl) delete obj;
          }

          cout << ANSI_DIM
               << "\nNOTE: Preselection fail counters are inclusive (one photon can increment multiple fail histograms).\n"
               << "Per-bin plots written under: " << outDir << "\n"
               << "3x3 table written to: " << JoinPath(outDir, "table3x3_preselectionFails.png") << "\n"
               << ANSI_RESET;
      }

      // =============================================================================
      // Section 3: general isolation QA
      // =============================================================================
      static void PrintIsoDecisionTable(Dataset& ds)
      {
        cout << ANSI_BOLD_CYN << "\n[SECTION 3] isoDecision table (" << ds.label << ")\n" << ANSI_RESET;

        const int wBin = 10;
        const int wN   = 14;

        cout << std::left
             << std::setw(wBin) << "pTbin"
             << std::right
             << std::setw(wN) << "N_iso"
             << std::setw(wN) << "N_nonIso"
             << std::setw(wN) << "isoFrac"
             << "\n";
        cout << string(wBin + 3*wN, '-') << "\n";

        const auto& bins = PtBins();
        for (int i = 0; i < kNPtBins; ++i)
        {
          const PtBin& b = bins[i];
          const string hname = "h_isoDecision" + b.suffix;

          TH1* h = GetObj<TH1>(ds, hname, true, true, false);
          double nIso = 0.0, nNon = 0.0;
          if (h)
          {
            nIso = h->GetBinContent(1);
            nNon = h->GetBinContent(2);
          }
          const double denom = nIso + nNon;
          const double frac  = (denom > 0.0) ? (nIso/denom) : 0.0;

          cout << std::left << std::setw(wBin) << b.label
               << std::right
               << std::setw(wN) << std::fixed << std::setprecision(0) << nIso
               << std::setw(wN) << nNon
               << std::setw(wN) << std::fixed << std::setprecision(4) << frac
               << "\n";
        }
      }

      // -------------------------------------------------------------------------
      // Helper: RECO / DATA isolation QA (3x3 tables + per-pT bin plots)
      //   - This is the "general isolation" part that runs for both data and sim.
      // -------------------------------------------------------------------------
      void RunIsolationQA_Reco(Dataset& ds, const string& outDir)
      {
        vector<string> common;
        common.push_back("Cuts: p_{T}^{#gamma} #geq 5 GeV, |#eta^{#gamma}| < 0.7");

        struct IsoHistDef { string base; string outStem; };
        const vector<IsoHistDef> isoHists = {
          {"h_Eiso",        "Eiso_total"},
          {"h_Eiso_emcal",  "Eiso_emcal"},
          {"h_Eiso_hcalin", "Eiso_hcalin"},
          {"h_Eiso_hcalout","Eiso_hcalout"}
        };

        for (const auto& def : isoHists)
        {
          // 3x3 table over pT bins
          Make3x3Table_TH1(ds, def.base, outDir,
                           string("table3x3_") + def.outStem + ".png",
                           "E_{iso} [GeV]", "Counts",
                           false, false, common);

          // Per pT bin plots
          for (int i = 0; i < kNPtBins; ++i)
          {
            const PtBin& b = PtBins()[i];
            const string hname = def.base + b.suffix;

            TH1* h = GetObj<TH1>(ds, hname, true, true, true);
            if (!h) continue;

            TH1* hc = CloneTH1(h, TString::Format("%s_%d", def.base.c_str(), i).Data());
            if (!hc) continue;

            const string fp = JoinPath(outDir, b.folder + "/" + def.outStem + "_" + b.folder + ".png");

            const string isoTitle = IsoTitleForBase(def.base);

            DrawAndSaveTH1_Iso(ds, hc, fp,
                                 "E_{iso} [GeV]", "Counts",
                                 isoTitle, b,
                                 false);

            delete hc;
          }
        }

          // -------------------------------------------------------------------------
          // (removed) ClusterIso comparison plot:
          //   h2_EisoBuilder_vs_ClusterIso<suffix>
          // PhotonClusterBuilder is now the ONLY isolation definition used.
          // -------------------------------------------------------------------------
      }

      // -------------------------------------------------------------------------
      // Helper: SIM-only TRUTH isolation QA (truth distributions + reco-inclusive
      // overlay + decision hist).
      // -------------------------------------------------------------------------
      void RunIsolationQA_TruthSim(Dataset& ds, const string& outDir)
      {
        const string truthDir = JoinPath(outDir, "TruthIsolation");
        EnsureDir(truthDir);

        // --- Read truth histograms (unsliced; live directly in /SIM/) ---
        TH1* hTruthIso = GetObj<TH1>(ds, "h_EisoTruth", true, true, true);
        TH1* hTruthDec = GetObj<TH1>(ds, "h_EisoTruthDecision", true, true, true);

        // --- Build reco-inclusive histogram by summing pT-sliced reco iso ---
        TH1* hRecoIsoSum = SumRecoPtSlicedTH1(ds, "h_Eiso");

        // 1) Truth isolation distribution (counts + shape)
        if (hTruthIso)
        {
          TH1* htCounts = CloneTH1(hTruthIso, "h_EisoTruth_counts_clone");
          TH1* htShape  = CloneTH1(hTruthIso, "h_EisoTruth_shape_clone");

          vector<string> lines = {
            "Truth isolation (prompt #gamma candidates in acceptance)",
            "Histogram: h_EisoTruth (SIM-only)"
          };

          if (htCounts)
          {
            DrawAndSaveTH1_Common(ds, htCounts,
              JoinPath(truthDir, "truthIso_counts.png"),
              "E_{T}^{iso,truth} [GeV]", "Counts", lines, false);
            delete htCounts;
          }

          if (htShape)
          {
            NormalizeToUnitArea(htShape);
            DrawAndSaveTH1_Common(ds, htShape,
              JoinPath(truthDir, "truthIso_shape.png"),
              "E_{T}^{iso,truth} [GeV]", "A.U.", lines, false);
            delete htShape;
          }
        }

        // 2) Truth isolation decision (PASS/FAIL)
        if (hTruthDec)
        {
          TH1* hd = CloneTH1(hTruthDec, "h_EisoTruthDecision_clone");
          vector<string> lines = {
            "Truth isolation decision",
            "Histogram: h_EisoTruthDecision (SIM-only)",
            "bin1=PASS, bin2=FAIL"
          };

          if (hd)
          {
            DrawAndSaveTH1_Common(ds, hd,
              JoinPath(truthDir, "truthIsoDecision.png"),
              "Decision", "Counts", lines, false);
            delete hd;
          }
        }

        // 3) Reco inclusive iso (counts + shape)
        if (hRecoIsoSum)
        {
          TH1* hrCounts = CloneTH1(hRecoIsoSum, "h_EisoRecoInclusive_counts_clone");
          TH1* hrShape  = CloneTH1(hRecoIsoSum, "h_EisoRecoInclusive_shape_clone");

          vector<string> lines = {
            "Reco isolation (inclusive over p_{T}^{#gamma} bins)",
            "Built offline by summing h_Eiso_pT_*"
          };

          if (hrCounts)
          {
            DrawAndSaveTH1_Common(ds, hrCounts,
              JoinPath(truthDir, "recoIsoInclusive_counts.png"),
              "E_{T}^{iso,reco} [GeV]", "Counts", lines, false);
            delete hrCounts;
          }

          if (hrShape)
          {
            NormalizeToUnitArea(hrShape);
            DrawAndSaveTH1_Common(ds, hrShape,
              JoinPath(truthDir, "recoIsoInclusive_shape.png"),
              "E_{T}^{iso,reco} [GeV]", "A.U.", lines, false);
            delete hrShape;
          }

          delete hRecoIsoSum;
          hRecoIsoSum = nullptr;
        }

        // 4) Overlay truth vs reco (shape) in common x-range
        // Choose overlap range automatically:
        //   reco typically ~[-5,12], truth often ~[0,50] => overlap ~[0,12].
        if (hTruthIso)
        {
          // Rebuild reco sum again for overlay (we deleted it above after saving).
          TH1* hRecoIsoSum2 = SumRecoPtSlicedTH1(ds, "h_Eiso");

          if (hRecoIsoSum2)
          {
            TH1* ht = CloneTH1(hTruthIso,   "h_EisoTruth_forOverlay");
            TH1* hr = CloneTH1(hRecoIsoSum2,"h_EisoReco_forOverlay");

            if (ht && hr)
            {
              // Determine overlap range from histogram axes
              const double xmin = std::max(ht->GetXaxis()->GetXmin(), hr->GetXaxis()->GetXmin());
              const double xmax = std::min(ht->GetXaxis()->GetXmax(), hr->GetXaxis()->GetXmax());

              // Guard against bad overlap
              double xlo = xmin;
              double xhi = xmax;
              if (!(std::isfinite(xlo) && std::isfinite(xhi) && xhi > xlo))
              {
                // fallback to a typical common window
                xlo = 0.0;
                xhi = 12.0;
              }

              // Normalize ONLY over the common visible range and draw the overlay
              NormalizeToUnitAreaInRange(ht, xlo, xhi);
              NormalizeToUnitAreaInRange(hr, xlo, xhi);

              // Draw only the overlap window for a fair visual comparison
              ht->GetXaxis()->SetRangeUser(xlo, xhi);
              hr->GetXaxis()->SetRangeUser(xlo, xhi);

              vector<string> extra = {
                "Overlay (shape): truth vs reco isolation",
                TString::Format("Normalized in-range: [%.2f, %.2f] GeV", xlo, xhi).Data(),
                "Truth: h_EisoTruth   |   Reco: sum of h_Eiso_pT_*"
              };

              DrawOverlayTwoTH1(ds, ht, hr,
                "Truth isolation", "Reco isolation (inclusive)",
                JoinPath(truthDir, "overlay_truthIso_vs_recoIso_shape.png"),
                "E_{T}^{iso} [GeV]", "A.U.", extra, false);
            }

            if (ht) delete ht;
            if (hr) delete hr;

            delete hRecoIsoSum2;
          }
        }

        // -------------------------------------------------------------------------
        // NEW (SIM only): reco isolation for TRUTH-ISOLATED signal photon matches
        //   Histogram per pT slice:
        //     h_EisoReco_truthSigMatched<suffix>
        //   Filled before tight/non-tight gating (SS-independent).
        // -------------------------------------------------------------------------
        {
          const string sigDir = JoinPath(truthDir, "TruthSigMatchedRecoIso");
          EnsureDir(sigDir);
          for (const auto& b : PtBins()) EnsureDir(JoinPath(sigDir, b.folder));

          for (int i = 0; i < kNPtBins; ++i)
          {
            const PtBin& b = PtBins()[i];
            const string hname = "h_EisoReco_truthSigMatched" + b.suffix;

            TH1* h = GetObj<TH1>(ds, hname, true, true, true);
            if (!h) continue;

            // Counts
            if (TH1* hc = CloneTH1(h, TString::Format("h_truthSigMatchedRecoIso_counts_%d", i).Data()))
            {
              const string fp = JoinPath(sigDir, b.folder + "/EisoReco_truthSigMatched_counts_" + b.folder + ".png");
              DrawAndSaveTH1_Iso(ds, hc, fp,
                                 "E_{T}^{iso,reco} [GeV]", "Counts",
                                 "E_{T}^{iso,reco} (truth iso signal match)", b,
                                 false);
              delete hc;
            }

            // Shape
            if (TH1* hs = CloneTH1(h, TString::Format("h_truthSigMatchedRecoIso_shape_%d", i).Data()))
            {
              NormalizeToUnitArea(hs);
              const string fp = JoinPath(sigDir, b.folder + "/EisoReco_truthSigMatched_shape_" + b.folder + ".png");
              DrawAndSaveTH1_Iso(ds, hs, fp,
                                 "E_{T}^{iso,reco} [GeV]", "A.U.",
                                 "E_{T}^{iso,reco} (truth iso signal match)", b,
                                 false);
              delete hs;
            }

            // Optional: overlay vs inclusive reco isolation in the same pT bin (shape)
            TH1* hReco = GetObj<TH1>(ds, "h_Eiso" + b.suffix, false, false, false);
            if (hReco)
            {
              TH1* hSigShape  = CloneTH1(h,     TString::Format("h_truthSigMatchedRecoIso_overlay_%d", i).Data());
              TH1* hRecoShape = CloneTH1(hReco, TString::Format("h_recoIso_overlay_%d", i).Data());

              if (hSigShape && hRecoShape)
              {
                NormalizeToUnitArea(hSigShape);
                NormalizeToUnitArea(hRecoShape);

                vector<string> lines = {
                  TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data(),
                  "Overlay (shape): truth-iso signal matched vs inclusive reco",
                  "Matched: h_EisoReco_truthSigMatched   |   Inclusive: h_Eiso"
                };

                const string fp = JoinPath(sigDir, b.folder + "/overlay_truthSigMatched_vs_inclusiveReco_shape_" + b.folder + ".png");

                DrawOverlayTwoTH1(ds, hSigShape, hRecoShape,
                  "Truth-iso signal matched", "Inclusive reco",
                  fp,
                  "E_{T}^{iso,reco} [GeV]", "A.U.", lines, false);
              }

              if (hSigShape)  delete hSigShape;
              if (hRecoShape) delete hRecoShape;
            }
          }
        }
      }

      void RunIsolationQA(Dataset& ds)
      {
        cout << ANSI_BOLD_CYN << "\n==============================\n"
             << "[SECTION 3] GENERAL ISOLATION QA (" << ds.label << ")\n"
             << "==============================" << ANSI_RESET << "\n";

        string outDir;
        if (ds.isSim) outDir = JoinPath(ds.outBase, "isoQAgeneral");
        else          outDir = JoinPath(ds.outBase, "baselineData/isoQAgeneral");

        EnsureDir(outDir);
        for (const auto& b : PtBins()) EnsureDir(JoinPath(outDir, b.folder));

        // 1) RECO/DATA isolation QA (always runs)
        RunIsolationQA_Reco(ds, outDir);

        // 2) Keep exactly where you had it before (after reco plots, before truth block)
        PrintIsoDecisionTable(ds);

        // ------------------------------------------------------------------
        // 2b) Photon-side ambiguity diagnostic:
        //     N = # of reco photon candidates that pass the SAME iso∧tight gate
        //         used for the leading photon selection in RecoilJets.cc.
        //     Filled once per event (only when a leading iso∧tight photon exists),
        //     and binned by the leading-photon pT bin.
        //
        // Hist family (written by RecoilJets.cc):
        //   h_nIsoTightPhoCand_pT_10_12, ..., h_nIsoTightPhoCand_pT_26_35
        // ------------------------------------------------------------------
        {
            const string outDirMult = JoinPath(outDir, "PhoCandMultiplicity");
            EnsureDir(outDirMult);

            vector<string> multLines = DefaultHeaderLines(ds);
            multLines.push_back("N_{#gamma}^{iso+tight} candidates per event");
            multLines.push_back("Counts ALL reco photon candidates passing the SAME iso#wedge tight gate");
            multLines.push_back("Filled once per event; pT bin = leading iso#wedge tight photon p_{T}");

            // 3x3 table of the N-candidate distribution in each pT bin
            Make3x3Table_TH1(ds,
                             "h_nIsoTightPhoCand",
                             outDirMult,
                             "table3x3_nIsoTightPhoCand.png",
                             "N_{#gamma}^{iso+tight} candidates",
                             "A.U.",
                             false,
                             true,
                             multLines);

            // Per-pT summary: fraction of events with N>=2
            const bool weighted = (ds.isSim && IsWeightedSIMSelected());

            cout << ANSI_BOLD_CYN
                 << "\n[PHOTON MULTIPLICITY] " << ds.label
                 << "  (hist base: h_nIsoTightPhoCand)\n"
                 << ANSI_RESET;

            const int wBin  = 10;
            const int wSum  = 16;
            const int wMean = 12;
            const int wFrac = 14;

            cout << std::left << std::setw(wBin) << "pTbin"
                 << std::right
                 << std::setw(wSum)  << (weighted ? "SumW" : "Nevents")
                 << std::setw(wMean) << "MeanN"
                 << std::setw(wFrac) << "Frac(N>=2)"
                 << "\n";
            cout << string(wBin + wSum + wMean + wFrac, '-') << "\n";

            vector<string> txt;
            txt.push_back("[PHOTON MULTIPLICITY] N iso+tight photon candidates per event");
            txt.push_back("Hist base: " + ds.topDirName + "/h_nIsoTightPhoCand_pT_*_*");
            txt.push_back("NOTE: Filled once per event when a leading iso+tight photon exists; pT bin is that leading photon.");
            if (weighted)
            {
              txt.push_back("NOTE: weighted merged SIM sample – SumW is sum of weights (not raw event count).");
            }
            txt.push_back("");
            txt.push_back(string("pTbin  ") + (weighted ? "SumW" : "Nevents") + "  MeanN  FracNge2");

            for (const auto& pb : PtBins())
            {
              const string hname = "h_nIsoTightPhoCand" + pb.suffix;
              TH1* h = GetObj<TH1>(ds, hname, true, false, false);

              double sumw   = 0.0;
              double meanN  = 0.0;
              double fracGe2= 0.0;

              if (h)
              {
                const int nb = h->GetNbinsX();
                sumw  = h->Integral(0, nb + 1);
                meanN = h->GetMean();

                const int b2 = h->GetXaxis()->FindBin(2.0); // N>=2
                const double ge2 = h->Integral(b2, nb + 1);
                fracGe2 = (sumw > 0.0) ? (ge2 / sumw) : 0.0;
              }

              std::ostringstream sSum;
              if (weighted) sSum << std::fixed << std::setprecision(6) << sumw;
              else          sSum << std::fixed << std::setprecision(0) << sumw;

              cout << std::left << std::setw(wBin) << pb.label
                   << std::right
                   << std::setw(wSum)  << sSum.str()
                   << std::setw(wMean) << std::fixed << std::setprecision(4) << meanN
                   << std::setw(wFrac) << std::fixed << std::setprecision(6) << fracGe2
                   << "\n";

              std::ostringstream line;
              line << pb.label << "  " << sSum.str()
                   << "  " << std::fixed << std::setprecision(4) << meanN
                   << "  " << std::fixed << std::setprecision(6) << fracGe2;
              txt.push_back(line.str());
            }

            WriteTextFile(JoinPath(outDirMult, "summary_nIsoTightPhoCand.txt"), txt);

            cout << ANSI_DIM
                 << "  -> Wrote: " << JoinPath(outDirMult, "table3x3_nIsoTightPhoCand.png") << "\n"
                 << "  -> Wrote: " << JoinPath(outDirMult, "summary_nIsoTightPhoCand.txt") << ANSI_RESET << "\n";
         }

         // 3) SIM-only truth isolation QA + overlays (only runs in SIM)
         if (ds.isSim)
         {
            RunIsolationQA_TruthSim(ds, outDir);
         }
      }

      // =============================================================================
      // Section 4: ABCD QA + purity + SS overlays + SIM leakage factors
      // =============================================================================
      void LoadLeakageFactorsFromSIM(Dataset& sim, LeakageFactors& lf)
      {
        lf.fB.assign(kNPtBins, 0.0);
        lf.fC.assign(kNPtBins, 0.0);
        lf.fD.assign(kNPtBins, 0.0);
        lf.available = false;

        if (!sim.isSim) return;

        bool any = false;

        cout << ANSI_BOLD_CYN << "\n[SECTION 4] Reading SIM leakage factors (h_sigABCD_MC_pT_*)\n" << ANSI_RESET;

        const int wBin = 10;
        const int wN   = 14;

        cout << std::left << std::setw(wBin) << "pTbin"
             << std::right
             << std::setw(wN) << "A_sig_MC"
             << std::setw(wN) << "B_sig_MC"
             << std::setw(wN) << "C_sig_MC"
             << std::setw(wN) << "D_sig_MC"
             << std::setw(wN) << "fB"
             << std::setw(wN) << "fC"
             << std::setw(wN) << "fD"
             << "\n";
        cout << string(wBin + 7*wN, '-') << "\n";

        for (int i = 0; i < kNPtBins; ++i)
        {
          const PtBin& b = PtBins()[i];
          const string hname = "h_sigABCD_MC" + b.suffix;

          TH1* h = GetObj<TH1>(sim, hname, true, true, false);
          double A = 0, B = 0, C = 0, D = 0;
          if (h)
          {
            A = h->GetBinContent(1);
            B = h->GetBinContent(2);
            C = h->GetBinContent(3);
            D = h->GetBinContent(4);
          }
          const double fB = (A > 0.0) ? (B/A) : 0.0;
          const double fC = (A > 0.0) ? (C/A) : 0.0;
          const double fD = (A > 0.0) ? (D/A) : 0.0;

          lf.fB[i] = fB;
          lf.fC[i] = fC;
          lf.fD[i] = fD;

          if (A > 0.0) any = true;

          cout << std::left << std::setw(wBin) << b.label
               << std::right
               << std::setw(wN) << std::fixed << std::setprecision(0) << A
               << std::setw(wN) << B
               << std::setw(wN) << C
               << std::setw(wN) << D
               << std::setw(wN) << std::fixed << std::setprecision(6) << fB
               << std::setw(wN) << fC
               << std::setw(wN) << fD
               << "\n";
        }

        lf.available = any;
        if (!lf.available)
        {
          cout << ANSI_BOLD_YEL << "[WARN] No nonzero A_sig_MC found; leakage correction will be unavailable.\n" << ANSI_RESET;
        }
      }
  
      static void Make3x3Table_ABCDCounts(Dataset& ds,
                                         const string& outDir,
                                         const string& outName,
                                         const vector<string>& commonLines)
      {
        EnsureDir(outDir);

        const char* xLabelsABCD[4] = {"N_{A}", "N_{B}", "N_{C}", "N_{D}"};
        const int   colorsABCD[4]  = {kGreen+2, kRed+1, kAzure+1, kMagenta+1};

        // Wider canvas for readability
        TCanvas c(
          TString::Format("c_abcd_cnt_tbl_%s", ds.label.c_str()).Data(),
          "c_abcd_cnt_tbl", 2400, 1400
        );
        c.Divide(3,2, 0.002, 0.002);

        vector<TObject*> keepAlive;
        keepAlive.reserve(kNPtBins * (1 + 4));

        const auto& bins = PtBins();

        for (int i = 0; i < kNPtBins; ++i)
        {
          c.cd(i+1);
          if (!gPad) continue;

          // Match your preselection-table style
          gPad->SetLeftMargin(0.17);
          gPad->SetRightMargin(0.08);
          gPad->SetTopMargin(0.14);
          gPad->SetBottomMargin(0.38);
          gPad->SetTicks(1,1);

          const PtBin& b   = bins[i];
          const string suf = b.suffix;

          // Read counts from the existing 1-bin histograms
          const double A = Read1BinCount(ds, "h_isIsolated_isTight"     + suf);
          const double B = Read1BinCount(ds, "h_notIsolated_isTight"    + suf);
          const double Cc = Read1BinCount(ds, "h_isIsolated_notTight"   + suf);
          const double D = Read1BinCount(ds, "h_notIsolated_notTight"   + suf);

          const double vals[4] = {A, B, Cc, D};

          double ymax = 0.0;
          for (int ib = 0; ib < 4; ++ib) ymax = std::max(ymax, vals[ib]);
          const double yMaxPlot = (ymax > 0.0) ? (1.35 * ymax) : 1.0;

          // Axis frame
          TH1F* hAxis = new TH1F(
            TString::Format("h_abcdCntAxis_%s_%s", ds.label.c_str(), b.folder.c_str()).Data(),
            "",
            4, 0.5, 4.5
          );
          hAxis->SetDirectory(nullptr);
          hAxis->SetStats(0);
          hAxis->SetMinimum(0.0);
          hAxis->SetMaximum(yMaxPlot);

          for (int ib = 1; ib <= 4; ++ib) hAxis->GetXaxis()->SetBinLabel(ib, xLabelsABCD[ib-1]);

          hAxis->GetYaxis()->SetTitle("Counts");
          hAxis->GetXaxis()->SetTitle("");
          hAxis->GetXaxis()->LabelsOption("h"); // horizontal labels; only 4 bins
          hAxis->GetXaxis()->SetLabelSize(0.080);
          hAxis->GetXaxis()->SetLabelOffset(0.012);

          hAxis->GetYaxis()->SetTitleSize(0.065);
          hAxis->GetYaxis()->SetTitleOffset(1.05);
          hAxis->GetYaxis()->SetLabelSize(0.055);

          hAxis->SetLineColor(1);
          hAxis->SetLineWidth(2);
          hAxis->SetFillStyle(0);
          hAxis->Draw("hist");

          // Bars (solid)
          for (int ib = 1; ib <= 4; ++ib)
          {
            TH1F* hb = new TH1F(
              TString::Format("h_abcdCntBar_%s_%s_b%d", ds.label.c_str(), b.folder.c_str(), ib).Data(),
              "",
              4, 0.5, 4.5
            );
            hb->SetDirectory(nullptr);
            hb->SetStats(0);

            hb->SetBinContent(ib, vals[ib-1]);
            hb->SetFillStyle(1001);
            hb->SetFillColor(colorsABCD[ib-1]);
            hb->SetLineColor(1);
            hb->SetLineWidth(2);

            hb->SetBarWidth(0.90);
            hb->SetBarOffset(0.05);

            hb->Draw("BAR SAME");
            keepAlive.push_back(hb);
          }

          // Numeric labels above bars
          TLatex t;
          t.SetTextFont(42);
          t.SetTextAlign(22);
          t.SetTextSize(0.060);

          for (int ib = 1; ib <= 4; ++ib)
          {
            const double y = vals[ib-1];
            if (y <= 0.0) continue;

            const double x = hAxis->GetXaxis()->GetBinCenter(ib);
            const double yText = std::min(y + 0.03*yMaxPlot, 0.95*yMaxPlot);
            t.DrawLatex(x, yText, TString::Format("%.0f", y).Data());
          }

          // Top annotations (NDC)
          TLatex tt;
          tt.SetTextFont(42);
          tt.SetNDC();

          // Title (top-right)
          tt.SetTextAlign(33);   // right, top
          tt.SetTextSize(0.062);
          tt.DrawLatex(0.95, 0.93, "ABCD counts");

          // pT bin (centered horizontally, same vertical band as title)
          tt.SetTextAlign(23);   // center, top
          tt.SetTextSize(0.070);
          tt.DrawLatex(0.50, 0.93,
              TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data()
          );

          keepAlive.push_back(hAxis);
        }

        SaveCanvas(c, JoinPath(outDir, outName));

        for (auto* obj : keepAlive) delete obj;
      }

      static void Make3x3Table_SSOverlay(Dataset& ds,
                                          const string& varKey,
                                          const string& outDir,
                                          const string& outName,
                                          const vector<string>& commonLines)
        {
          EnsureDir(outDir);

          const vector<string> regionTags = {
            "isIsolated_isTight",
            "notIsolated_isTight",
            "isIsolated_notTight",
            "notIsolated_notTight"
          };

          // Use a unique canvas name to avoid ROOT name collisions across calls
          TCanvas c(
            TString::Format("c_ss_tbl_%s_%s", ds.label.c_str(), varKey.c_str()).Data(),
            "c_ss_tbl", 1500, 1200
          );
          c.Divide(3,2, 0.001, 0.001);

          // Keep cloned histograms alive until AFTER SaveCanvas()
          vector<TH1*> keepAlive;
          keepAlive.reserve(kNPtBins * 4);

          const auto& bins = PtBins();

          for (int i = 0; i < kNPtBins; ++i)
          {
            c.cd(i+1);
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.10);

            const PtBin& b = bins[i];

            vector<TH1*> hc(4, nullptr);

            for (int r = 0; r < 4; ++r)
            {
              const string hname = "h_ss_" + varKey + "_" + regionTags[r] + b.suffix;

              TH1* h = GetObj<TH1>(ds, hname, true, true, true);
              if (!h) continue;

              hc[r] = CloneTH1(h, TString::Format("ss_%s_%d_%d", varKey.c_str(), i, r).Data());
              if (!hc[r]) continue;

              NormalizeToUnitArea(hc[r]);
              hc[r]->SetLineWidth(2);
              hc[r]->SetLineColor((r==0)?1:(r==1)?2:(r==2)?4:6);

              // IMPORTANT: do NOT delete until after SaveCanvas
              keepAlive.push_back(hc[r]);
            }

            TH1* first = nullptr;
            for (int r = 0; r < 4; ++r) if (hc[r]) { first = hc[r]; break; }

            if (!first)
            {
              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextSize(0.06);
              t.DrawLatex(0.15, 0.55, "MISSING");
              t.SetTextSize(0.05);
              std::ostringstream s;
              s << "SS " << varKey << "  pT^{#gamma}: " << b.lo << "-" << b.hi;
              t.DrawLatex(0.15, 0.45, s.str().c_str());
              continue;
            }

            double ymax = 0.0;
            for (auto* p : hc) if (p) ymax = std::max(ymax, p->GetMaximum());
            first->SetMaximum(ymax * 1.25);

            first->SetTitle("");
            first->GetXaxis()->SetTitle(varKey.c_str());
            first->GetYaxis()->SetTitle("A.U.");
            first->Draw("hist");
            for (int r = 0; r < 4; ++r) if (hc[r] && hc[r] != first) hc[r]->Draw("hist same");

            vector<string> lines = commonLines;
            lines.push_back(TString::Format("SS %s  pT^{#gamma}: %d-%d GeV", varKey.c_str(), b.lo, b.hi).Data());
            DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);
          }

          // Now the canvas can repaint with valid objects
          SaveCanvas(c, JoinPath(outDir, outName));

          // Cleanup AFTER saving
          for (auto* h : keepAlive) delete h;
      }

      void RunABCDPurityAndSidebandSubtraction(Dataset& ds, const LeakageFactors& lf)
      {
        cout << ANSI_BOLD_CYN << "\n==============================\n"
             << "[SECTION 4] ABCD QA + PURITY + SS (" << ds.label << ")\n"
             << "==============================" << ANSI_RESET << "\n";

        string outDir;
        if (ds.isSim) outDir = JoinPath(ds.outBase, "PurityABCD");
        else          outDir = JoinPath(ds.outBase, "baselineData/PurityABCD");
        EnsureDir(outDir);

        vector<double> purityRaw(kNPtBins, 0.0);
        vector<double> purityCorr(kNPtBins, 0.0);
        vector<bool>   hasCorr(kNPtBins, false);

        cout << ANSI_BOLD_CYN << "\n[ABCD COUNTS + PURITY] " << ds.label << "\n" << ANSI_RESET;

        const int wBin = 10;
        const int wN   = 12;

        cout << std::left << std::setw(wBin) << "pTbin"
             << std::right
             << std::setw(wN) << "A"
             << std::setw(wN) << "B"
             << std::setw(wN) << "C"
             << std::setw(wN) << "D"
             << std::setw(wN) << "A_sig"
             << std::setw(wN) << "Pur_raw"
             << std::setw(wN) << "Pur_corr"
             << "\n";
        cout << string(wBin + 7*wN, '-') << "\n";

        for (int i = 0; i < kNPtBins; ++i)
        {
            const PtBin& b = PtBins()[i];
            const string& suf = b.suffix;

            const double A = Read1BinCount(ds, "h_isIsolated_isTight" + suf);
            const double B = Read1BinCount(ds, "h_notIsolated_isTight" + suf);
            const double C = Read1BinCount(ds, "h_isIsolated_notTight" + suf);
            const double D = Read1BinCount(ds, "h_notIsolated_notTight" + suf);

            double Asig = 0.0;
            double Praw = 0.0;
            if (A > 0.0 && D > 0.0)
            {
              Asig = A - B*(C/D);
              if (Asig < 0.0) Asig = 0.0;
              Praw = Asig / A;
            }

            purityRaw[i] = Praw;

            double Pc = Praw;
            bool okCorr = false;
            if (lf.available)
            {
              double SA = 0.0;
              const bool ok = SolveLeakageCorrectedSA(A,B,C,D, lf.fB[i], lf.fC[i], lf.fD[i], SA);
              if (ok && A > 0.0)
              {
                Pc = SA / A;
                okCorr = true;
              }
              else
              {
                Pc = Praw;
                okCorr = false;
              }
            }

            purityCorr[i] = Pc;
            hasCorr[i]    = okCorr;

            cout << std::left << std::setw(wBin) << b.label
                 << std::right
                 << std::setw(wN) << std::fixed << std::setprecision(0) << A
                 << std::setw(wN) << B
                 << std::setw(wN) << C
                 << std::setw(wN) << D
                 << std::setw(wN) << Asig
                 << std::setw(wN) << std::fixed << std::setprecision(4) << Praw
                 << std::setw(wN) << std::fixed << std::setprecision(4) << Pc
                 << "\n";

            if (lf.available && !okCorr)
            {
              cout << ANSI_BOLD_YEL << "  [WARN] leakage solver failed for pT " << b.label
                   << " (fallback to raw)" << ANSI_RESET << "\n";
            }
          }

          // ---------------------------------------------------------------------------
          // 3x3 table of ABCD region counts (N_A, N_B, N_C, N_D) vs pT bin
          // ---------------------------------------------------------------------------
          {
            const string cntDir = JoinPath(outDir, "ABCDCounts");
            EnsureDir(cntDir);

            vector<string> cntCommon;
            cntCommon.push_back("ABCD regions (preselection pass)");
            cntCommon.push_back("A=iso&tight   B=nonIso&tight");
            cntCommon.push_back("C=iso&nonTight   D=nonIso&nonTight");

            Make3x3Table_ABCDCounts(ds, cntDir, "table3x3_ABCDCounts.png", cntCommon);
          }

          // purity_raw plot
          {
            TH1F hPur("hPurRaw","hPurRaw", kNPtBins, 0.5, kNPtBins + 0.5);
          hPur.SetDirectory(nullptr);
          hPur.GetYaxis()->SetRangeUser(0.0, 1.05);
          hPur.GetXaxis()->SetTitle("p_{T}^{#gamma} bin");
          hPur.GetYaxis()->SetTitle("Purity (raw ABCD)");
          for (int i = 0; i < kNPtBins; ++i)
          {
            hPur.SetBinContent(i+1, purityRaw[i]);
            hPur.GetXaxis()->SetBinLabel(i+1, PtBins()[i].label.c_str());
          }

          TCanvas c("c_pur_raw","c_pur_raw",900,700);
          ApplyCanvasMargins1D(c);

          hPur.SetLineWidth(2);
          hPur.SetMarkerStyle(20);
          hPur.SetMarkerSize(1.2);
          hPur.Draw("E1");

          vector<string> box;
          box.push_back("ABCD purity (raw)");
          box.push_back("Cuts: p_{T}^{#gamma} #geq 5 GeV, |#eta| < 0.7, preselection pass");
          box.push_back("Iso: E_{iso} < 1.08128 + 0.0299107 p_{T}^{#gamma}");
          box.push_back("NonIso: E_{iso} > isoThresh + 1 GeV");
          DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
          DrawLatexLines(0.14, 0.78, box, 0.030, 0.040);

          const string fp = JoinPath(outDir, ds.isSim ? "purity_raw_SIM.png" : "purity_raw_DATA.png");
          SaveCanvas(c, fp);
        }

        // overlay raw vs corrected if correction used
        bool anyCorr = false;
        for (bool b : hasCorr) if (b) { anyCorr = true; break; }

        if (anyCorr)
        {
          TH1F hR("hPurRaw2","hPurRaw2", kNPtBins, 0.5, kNPtBins + 0.5);
          TH1F hC("hPurCor2","hPurCor2", kNPtBins, 0.5, kNPtBins + 0.5);
          hR.SetDirectory(nullptr);
          hC.SetDirectory(nullptr);

          hR.GetYaxis()->SetRangeUser(0.0, 1.05);
          hR.GetXaxis()->SetTitle("p_{T}^{#gamma} bin");
          hR.GetYaxis()->SetTitle("Purity");

          for (int i = 0; i < kNPtBins; ++i)
          {
            hR.SetBinContent(i+1, purityRaw[i]);
            hC.SetBinContent(i+1, purityCorr[i]);
            hR.GetXaxis()->SetBinLabel(i+1, PtBins()[i].label.c_str());
            hC.GetXaxis()->SetBinLabel(i+1, PtBins()[i].label.c_str());
          }

          hR.SetLineWidth(2); hC.SetLineWidth(2);
          hR.SetMarkerStyle(20); hC.SetMarkerStyle(24);
          hR.SetLineColor(1); hC.SetLineColor(2);
          hR.SetMarkerColor(1); hC.SetMarkerColor(2);

          TCanvas c("c_pur_ov","c_pur_ov",900,700);
          ApplyCanvasMargins1D(c);

          hR.Draw("E1");
          hC.Draw("E1 same");

          TLegend leg(0.62,0.77,0.92,0.90);
          leg.SetTextFont(42);
          leg.SetTextSize(0.033);
          leg.AddEntry(&hR, "Raw ABCD", "lp");
          leg.AddEntry(&hC, "Leakage-corrected", "lp");
          leg.Draw();

          vector<string> box;
          box.push_back("ABCD purity: raw vs leakage-corrected");
          if (!ds.isSim) box.push_back("Leakage factors from SIM h_{sigABCD}^{MC}");
          DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
          DrawLatexLines(0.14, 0.78, box, 0.030, 0.040);

          const string fp = JoinPath(outDir, ds.isSim
            ? "purity_raw_vs_leakageCorrected_SIM.png"
            : "purity_raw_vs_leakageCorrected_DATA.png");

          SaveCanvas(c, fp);
        }

        // SS overlays
        cout << ANSI_BOLD_CYN << "\n[SS OVERLAYS] " << ds.label << "\n" << ANSI_RESET;

        const string ssBase = JoinPath(outDir, "SS");
        EnsureDir(ssBase);

        const vector<string> varKeys = {"weta","wphi","e11e33","e32e35","et1"};
        vector<string> ssCommon;
        ssCommon.push_back("SS overlays: preselection pass");
        ssCommon.push_back("Tight: 0 fails; NonTight: #geq2 fails; (1-fail excluded)");
        ssCommon.push_back("Iso: E_{iso} < isoThresh; NonIso: E_{iso} > isoThresh + 1 GeV");

        for (const auto& var : varKeys)
        {
          const string varDir = JoinPath(ssBase, var);
          EnsureDir(varDir);

          Make3x3Table_SSOverlay(ds, var, varDir, string("table3x3_ss_") + var + ".png", ssCommon);

          // per pT-bin overlay
          for (int i = 0; i < kNPtBins; ++i)
          {
            const PtBin& b = PtBins()[i];
            const string pDir = JoinPath(varDir, b.folder);
            EnsureDir(pDir);

            const vector<string> regionTags = {
              "isIsolated_isTight",
              "notIsolated_isTight",
              "isIsolated_notTight",
              "notIsolated_notTight"
            };
            const vector<string> regionLabels = {
              "A: iso&tight",
              "B: nonIso&tight",
              "C: iso&nonTight",
              "D: nonIso&nonTight"
            };
            const int colors[4] = {1,2,4,6};

            vector<TH1*> hcl(4, nullptr);
            for (int r = 0; r < 4; ++r)
            {
              const string hname = "h_ss_" + var + "_" + regionTags[r] + b.suffix;
              TH1* h = GetObj<TH1>(ds, hname, true, true, true);
              if (!h) continue;
              hcl[r] = CloneTH1(h, TString::Format("ss_%s_%d_%d", var.c_str(), i, r).Data());
              if (hcl[r])
              {
                NormalizeToUnitArea(hcl[r]);
                hcl[r]->SetLineWidth(2);
                hcl[r]->SetLineColor(colors[r]);
              }
            }

            TH1* first = nullptr;
            for (int r = 0; r < 4; ++r) if (hcl[r]) { first = hcl[r]; break; }
            if (!first)
            {
              for (auto* p : hcl) if (p) delete p;
              continue;
            }

            double ymax = 0.0;
            for (auto* p : hcl) if (p) ymax = std::max(ymax, p->GetMaximum());

            TCanvas c("c_ss","c_ss",900,700);
            ApplyCanvasMargins1D(c);

            first->SetTitle("");
            first->GetXaxis()->SetTitle(var.c_str());
            first->GetYaxis()->SetTitle("A.U.");
            first->SetMaximum(ymax * 1.25);

            first->Draw("hist");
            for (int r = 0; r < 4; ++r) if (hcl[r] && hcl[r] != first) hcl[r]->Draw("hist same");

            TLegend leg(0.58,0.70,0.92,0.90);
            leg.SetTextFont(42);
            leg.SetTextSize(0.030);
            for (int r = 0; r < 4; ++r) if (hcl[r]) leg.AddEntry(hcl[r], regionLabels[r].c_str(), "l");
            leg.Draw();

            DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);

            vector<string> box = ssCommon;
            box.insert(box.begin(), TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
            DrawLatexLines(0.14,0.78, box, 0.028, 0.038);

            SaveCanvas(c, JoinPath(pDir, string("ss_") + var + "_" + b.folder + ".png"));

            for (auto* p : hcl) if (p) delete p;
          }
        }
      }

      // =============================================================================
      // Section 5: Jet QA + recoil jet QA suite
      // =============================================================================
      static void PlotJetQA_AllOrIncl(Dataset& ds,
                                     const string& baseOut,
                                     const string& modeTag, // "all" or "incl"
                                     bool fiducialLinesForEta,
                                     bool includeEventLevel)
      {
        EnsureDir(baseOut);

        const double nEvtFallback = ReadEventCount(ds);

        for (const auto& rKey : kRKeys)
        {
          const double R = RFromKey(rKey);
          const double etaFidAbs = FidEtaAbsFromKey(rKey);

          const string rOut = JoinPath(baseOut, rKey);
          const string dirShape    = JoinPath(rOut, "shape");
          const string dirJpe      = JoinPath(rOut, "jetsPerEvent");
          const string dir2DJpe    = JoinPath(rOut, "2D_jetsPerEvent");
          const string dirProfiles = JoinPath(rOut, "profiles");
          const string dirEvent    = JoinPath(rOut, "event");

          EnsureDir(dirShape);
          EnsureDir(dirJpe);
          EnsureDir(dir2DJpe);
          EnsureDir(dirProfiles);
          if (includeEventLevel) EnsureDir(dirEvent);

          const double Nevt = DetermineNevtForRKey(ds, rKey, nEvtFallback);

          vector<string> notes;
          if (Nevt <= 0.0) notes.push_back("Nevt is 0; jetsPerEvent scaling is meaningless/skipped.");

          auto H1 = [&](const string& stem)->TH1* {
            return GetObj<TH1>(ds, stem + "_" + modeTag + "_" + rKey, true, true, true);
          };
          auto H2 = [&](const string& stem)->TH2* {
            return GetObj<TH2>(ds, stem + "_" + modeTag + "_" + rKey, true, true, true);
          };

          TH1* hPt       = H1("h_jetPt");
          TH1* hEta      = H1("h_jetEta");
          TH1* hPhi      = H1("h_jetPhi");
          TH2* hEtaPhi   = H2("h_jetEtaPhi");
          TH1* hMass     = H1("h_jetMass");
          TH2* hMassVsPt = H2("h_jetMassVsPt");

          auto doShapeAndJpe = [&](TH1* hIn, const string& stemOut, const string& xTitle,
                                   bool logy, bool drawFid = false)
          {
            if (!hIn) { notes.push_back("Missing " + stemOut + " (" + modeTag + ")"); return; }

            TH1* hShape = CloneTH1(hIn, stemOut + "_shape");
            TH1* hJpe   = CloneTH1(hIn, stemOut + "_jpe");
            if (!hShape || !hJpe) { if (hShape) delete hShape; if (hJpe) delete hJpe; return; }

            NormalizeToUnitArea(hShape);

            vector<string> lines;
            lines.push_back(string("Jet QA: ") + modeTag + "  " + rKey + TString::Format(" (R=%.1f)", R).Data());
            lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
            if (modeTag == "incl") lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());

            DrawAndSaveTH1_Common(ds, hShape,
              JoinPath(dirShape, stemOut + "_shape.png"),
              xTitle, "A.U.", lines, logy, drawFid, etaFidAbs);

            if (Nevt > 0.0) hJpe->Scale(1.0/Nevt);
            vector<string> lines2 = lines;
            lines2.push_back(TString::Format("Scaled: 1/N_{evt} (N_{evt}=%.0f)", Nevt).Data());

            DrawAndSaveTH1_Common(ds, hJpe,
              JoinPath(dirJpe, stemOut + "_jetsPerEvent.png"),
              xTitle, "Jets / event / bin", lines2, logy, drawFid, etaFidAbs);

            delete hShape;
            delete hJpe;
          };

          doShapeAndJpe(hPt,   "jetPt",   "p_{T}^{jet} [GeV]", true,  false);
          doShapeAndJpe(hEta,  "jetEta",  "#eta_{jet}",        false, fiducialLinesForEta);
          doShapeAndJpe(hPhi,  "jetPhi",  "#phi_{jet}",        false, false);
          doShapeAndJpe(hMass, "jetMass", "m_{jet} [GeV]",     false, false);

          // eta-phi occupancy (jets/event only)
          if (hEtaPhi)
          {
            TH2* h2 = CloneTH2(hEtaPhi, "jetEtaPhi_jpe");
            if (h2)
            {
              if (Nevt > 0.0) h2->Scale(1.0/Nevt);

              vector<string> lines;
              lines.push_back(string("Jet #eta-#phi occupancy: ") + modeTag + "  " + rKey);
              lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
              if (modeTag == "incl") lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());
              lines.push_back(TString::Format("Scaled: 1/N_{evt} (N_{evt}=%.0f)", Nevt).Data());

              DrawAndSaveTH2_Common(ds, h2,
                JoinPath(dir2DJpe, "jetEtaPhi_jetsPerEvent.png"),
                "#eta_{jet}", "#phi_{jet}", "Jets / event / bin", lines,
                false, fiducialLinesForEta, etaFidAbs);

              delete h2;
            }
          }
          else notes.push_back("Missing h_jetEtaPhi_" + modeTag + "_" + rKey);

          // mass vs pt (jets/event) + ProfileX
          if (hMassVsPt)
          {
            TH2* h2 = CloneTH2(hMassVsPt, "jetMassVsPt_jpe");
            if (h2)
            {
              if (Nevt > 0.0) h2->Scale(1.0/Nevt);

              vector<string> lines;
              lines.push_back(string("Jet mass vs p_{T}: ") + modeTag + "  " + rKey);
              lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
              if (modeTag == "incl") lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());
              lines.push_back(TString::Format("Scaled: 1/N_{evt} (N_{evt}=%.0f)", Nevt).Data());

              DrawAndSaveTH2_Common(ds, h2,
                JoinPath(dir2DJpe, "jetMassVsPt_jetsPerEvent.png"),
                "p_{T}^{jet} [GeV]", "m_{jet} [GeV]", "Jets / event / bin",
                lines, false);

              TProfile* p = h2->ProfileX("p_mass_vs_pt");
              if (p)
              {
                p->SetDirectory(nullptr);
                TH1* asH = dynamic_cast<TH1*>(p);
                vector<string> l2 = lines;
                l2.push_back("ProfileX: mean m_{jet} vs p_{T}");
                DrawAndSaveTH1_Common(ds, asH,
                  JoinPath(dirProfiles, "profile_meanJetMass_vs_jetPt.png"),
                  "p_{T}^{jet} [GeV]", "<m_{jet}> [GeV]", l2, false);
                delete p;
              }
              delete h2;
            }
          }
          else notes.push_back("Missing h_jetMassVsPt_" + modeTag + "_" + rKey);

          // event-level (incl only)
          map<string,double> scalars;
          scalars["R"] = R;
          scalars["etaFidAbs"] = etaFidAbs;
          scalars["Nevt"] = Nevt;

          if (includeEventLevel)
          {
            auto H1evt = [&](const string& stem)->TH1* {
              return GetObj<TH1>(ds, stem + "_" + rKey, true, true, true);
            };

            TH1* hNJets   = H1evt("h_nJets");
            TH1* hHT      = H1evt("h_HT");
            TH1* hLeadPt  = H1evt("h_leadJetPt");
            TH1* hLeadEta = H1evt("h_leadJetEta");
            TH1* hLeadPhi = H1evt("h_leadJetPhi");
            TH1* hSubPt   = H1evt("h_subleadJetPt");

            auto doEventShape = [&](TH1* hIn, const string& stemOut, const string& xTitle,
                                    bool logy, bool drawFid = false)
            {
              if (!hIn) { notes.push_back("Missing " + stemOut + "_" + rKey); return; }
              TH1* hc = CloneTH1(hIn, stemOut + "_shape");
              if (!hc) return;
              NormalizeToUnitArea(hc);

              vector<string> lines;
              lines.push_back(string("Event-level (fid jets): ") + rKey);
              lines.push_back("Jets: p_{T}^{jet} #geq 10 GeV");
              lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());

              DrawAndSaveTH1_Common(ds, hc,
                JoinPath(dirEvent, stemOut + "_shape.png"),
                xTitle, "A.U.", lines, logy, drawFid, etaFidAbs);

              delete hc;
            };

            doEventShape(hNJets,   "nJets",       "N_{jets}^{fid}",              false, false);
            doEventShape(hHT,      "HT",          "H_{T} [GeV]",                 false, false);
            doEventShape(hLeadPt,  "leadJetPt",   "p_{T}^{lead} [GeV]",          true,  false);
            doEventShape(hLeadEta, "leadJetEta",  "#eta_{lead}",                 false, true);
            doEventShape(hLeadPhi, "leadJetPhi",  "#phi_{lead}",                 false, false);
            doEventShape(hSubPt,   "subleadJetPt","p_{T}^{sublead} [GeV]",       true,  false);

            const double f_ge1 = (Nevt > 0.0 && hLeadPt) ? (hLeadPt->GetEntries()/Nevt) : 0.0;
            const double f_ge2 = (Nevt > 0.0 && hSubPt)  ? (hSubPt->GetEntries()/Nevt)  : 0.0;
            scalars["f_ge1"] = f_ge1;
            scalars["f_ge2"] = f_ge2;
          }

          // jets/event summary scalars
          if (hPt)
          {
            const double nJets = hPt->Integral(0, hPt->GetNbinsX()+1);
            scalars["Njets_total"] = nJets;
            scalars["jets_per_event"] = (Nevt > 0.0) ? (nJets / Nevt) : 0.0;
            scalars["mean_jetPt"] = hPt->GetMean();
          }
          if (hMass) scalars["mean_jetMass"] = hMass->GetMean();

          WriteJetSummaryTxt(JoinPath(rOut, "summary.txt"), rKey, Nevt, scalars, notes);
        }

        // r02 vs r04 overlays for jetPt/eta/phi shapes
        {
          const string overDir = JoinPath(baseOut, "overlays");
          EnsureDir(overDir);

          auto getH = [&](const string& stem, const string& mode, const string& rKey)->TH1*
          {
            return GetObj<TH1>(ds, stem + "_" + mode + "_" + rKey, true, true, true);
          };

          TH1* pt02  = getH("h_jetPt", modeTag, "r02");
          TH1* pt04  = getH("h_jetPt", modeTag, "r04");
          TH1* eta02 = getH("h_jetEta", modeTag, "r02");
          TH1* eta04 = getH("h_jetEta", modeTag, "r04");
          TH1* phi02 = getH("h_jetPhi", modeTag, "r02");
          TH1* phi04 = getH("h_jetPhi", modeTag, "r04");

          if (pt02 && pt04)
          {
            TH1* a = CloneTH1(pt02, "pt02_shape"); TH1* b = CloneTH1(pt04, "pt04_shape");
            NormalizeToUnitArea(a); NormalizeToUnitArea(b);
            DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
              JoinPath(overDir, "overlay_jetPt_" + modeTag + "_shape_logy.png"),
              "p_{T}^{jet} [GeV]", "A.U.",
              {string("Overlay: ") + modeTag + " jet p_{T} (shape)"},
              true);
            delete a; delete b;
          }
          if (eta02 && eta04)
          {
            TH1* a = CloneTH1(eta02, "eta02_shape"); TH1* b = CloneTH1(eta04, "eta04_shape");
            NormalizeToUnitArea(a); NormalizeToUnitArea(b);
            DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
              JoinPath(overDir, "overlay_jetEta_" + modeTag + "_shape.png"),
              "#eta_{jet}", "A.U.",
              {string("Overlay: ") + modeTag + " jet #eta (shape)"},
              false);
            delete a; delete b;
          }
          if (phi02 && phi04)
          {
            TH1* a = CloneTH1(phi02, "phi02_shape"); TH1* b = CloneTH1(phi04, "phi04_shape");
            NormalizeToUnitArea(a); NormalizeToUnitArea(b);
            DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
              JoinPath(overDir, "overlay_jetPhi_" + modeTag + "_shape.png"),
              "#phi_{jet}", "A.U.",
              {string("Overlay: ") + modeTag + " jet #phi (shape)"},
              false);
            delete a; delete b;
          }
        }
      }

      void RunGeneralJetQA(Dataset& ds)
      {
        cout << ANSI_BOLD_CYN << "\n==============================\n"
             << "[SECTION 5A/5B] GeneralJetQA (" << ds.label << ")\n"
             << "==============================" << ANSI_RESET << "\n";

        const string baseOut = ds.isSim ? ds.outBase : JoinPath(ds.outBase, "baselineData");

        // 5A: pTcut_noFiducial (all)
        PlotJetQA_AllOrIncl(ds, JoinPath(baseOut, "GeneralJetQA/pTcut_noFiducial"),
                             "all", true, false);

        // 5B: fiducial inclusive jets (incl) + event-level
        const string baseIncl = JoinPath(baseOut, "GeneralJetQA/pTcutFiducialJets");
        PlotJetQA_AllOrIncl(ds, baseIncl, "incl", true, true);

        // extra overlays for incl: jetPt_incl, nJets, HT
        {
          const string overDir = JoinPath(baseIncl, "overlays");
          EnsureDir(overDir);

          TH1* n02  = GetObj<TH1>(ds, "h_nJets_r02", true, true, true);
          TH1* n04  = GetObj<TH1>(ds, "h_nJets_r04", true, true, true);
          TH1* ht02 = GetObj<TH1>(ds, "h_HT_r02", true, true, true);
          TH1* ht04 = GetObj<TH1>(ds, "h_HT_r04", true, true, true);
          TH1* pt02 = GetObj<TH1>(ds, "h_jetPt_incl_r02", true, true, true);
          TH1* pt04 = GetObj<TH1>(ds, "h_jetPt_incl_r04", true, true, true);

          if (pt02 && pt04)
          {
            TH1* a = CloneTH1(pt02, "pt_incl_r02_shape"); TH1* b = CloneTH1(pt04, "pt_incl_r04_shape");
            NormalizeToUnitArea(a); NormalizeToUnitArea(b);
            DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
              JoinPath(overDir, "overlay_jetPt_incl_shape_logy.png"),
              "p_{T}^{jet} [GeV]", "A.U.",
              {"Overlay: fiducial inclusive jet p_{T} (shape)"},
              true);
            delete a; delete b;
          }
          if (n02 && n04)
          {
            TH1* a = CloneTH1(n02, "nJets_r02_shape"); TH1* b = CloneTH1(n04, "nJets_r04_shape");
            NormalizeToUnitArea(a); NormalizeToUnitArea(b);
            DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
              JoinPath(overDir, "overlay_nJets_shape.png"),
              "N_{jets}^{fid}", "A.U.",
              {"Overlay: N_{jets}^{fid} (shape)"},
              false);
            delete a; delete b;
          }
          if (ht02 && ht04)
          {
            TH1* a = CloneTH1(ht02, "HT_r02_shape"); TH1* b = CloneTH1(ht04, "HT_r04_shape");
            NormalizeToUnitArea(a); NormalizeToUnitArea(b);
            DrawOverlayTwoTH1(ds, a, b, "r02 (R=0.2)", "r04 (R=0.4)",
              JoinPath(overDir, "overlay_HT_shape.png"),
              "H_{T} [GeV]", "A.U.",
              {"Overlay: H_{T} (shape)"},
              false);
            delete a; delete b;
          }
        }
      }

      void RunMatchQA(Dataset& ds, MatchCache& mc)
      {
            cout << ANSI_BOLD_CYN << "\n==============================\n"
                 << "[SECTION 5C] #gamma-jet MatchQA (" << ds.label << ")\n"
                 << "==============================" << ANSI_RESET << "\n";

            // Reset/prepare cache (used later for r02/r04 overlays)
            InitMatchCache(mc);

            // -------------------------------------------------------------------------
            // Output base directory
            // -------------------------------------------------------------------------
            string baseOut;
            if (ds.isSim) baseOut = JoinPath(ds.outBase, "RecoilJetQA/MatchQA");
            else          baseOut = JoinPath(ds.outBase, "baselineData/RecoilJetQA/MatchQA");
            EnsureDir(baseOut);

            // -------------------------------------------------------------------------
            // Small structs for cleanliness
            // -------------------------------------------------------------------------
            struct MatchDirs
            {
              string rOut;
              string dir2D;
              string dirProj;
            };

            auto MakeDirsForRKey = [&](const string& rKey)->MatchDirs
            {
              MatchDirs D;
              D.rOut    = JoinPath(baseOut, rKey);
              D.dir2D   = JoinPath(D.rOut, "2D");
              D.dirProj = JoinPath(D.rOut, "projections");
              EnsureDir(D.rOut);
              EnsureDir(D.dir2D);
              EnsureDir(D.dirProj);
              return D;
            };

            auto RLabel = [&](const string& rKey)->string
            {
              return TString::Format(" (R=%.1f)", RFromKey(rKey)).Data();
            };

            // -------------------------------------------------------------------------
            // Helper: draw match status TH2 + compute per-pT fractions + fill MatchCache
            // -------------------------------------------------------------------------
            auto HandleMatchStatus =
              [&](const string& rKey,
                  const MatchDirs& D,
                  vector<string>& notes)->TH2*
            {
              TH2* hStatus = GetObj<TH2>(ds, "h_match_status_vs_pTgamma_" + rKey, true, true, true);
              if (!hStatus)
              {
                notes.push_back("Missing h_match_status_vs_pTgamma_" + rKey);
                return nullptr;
              }

              // 2D plot (existing behavior)
              {
                TH2* hc = CloneTH2(hStatus, "status_clone");
                DrawAndSaveTH2_Common(ds, hc,
                  JoinPath(D.dir2D, "match_status_vs_pTgamma.png"),
                  "p_{T}^{#gamma} [GeV]", "Status", "Counts",
                  {string("Match status vs p_{T}^{#gamma}"), rKey + RLabel(rKey)},
                  false);
                delete hc;
              }

              // Compute per-bin fractions from status bins
              // Status binning (Y):
              //   1 = NoJetPt
              //   2 = NoJetEta
              //   3 = NoBackToBack
              //   4 = Matched
              vector<double> x(kNPtBins, 0.0);
              vector<double> f1(kNPtBins, 0.0), f2(kNPtBins, 0.0), f3(kNPtBins, 0.0), f4(kNPtBins, 0.0);

              double totAll = 0.0, tot1=0.0, tot2=0.0, tot3=0.0, tot4=0.0;

              for (int i = 0; i < kNPtBins; ++i)
              {
                const PtBin& b = PtBins()[i];
                x[i] = 0.5*(b.lo + b.hi);

                const int xbin = i+1;
                const double n1 = hStatus->GetBinContent(xbin, 1);
                const double n2 = hStatus->GetBinContent(xbin, 2);
                const double n3 = hStatus->GetBinContent(xbin, 3);
                const double n4 = hStatus->GetBinContent(xbin, 4);

                const double nLead = n1+n2+n3+n4;
                mc.NphoLead[rKey][i]    = nLead;
                mc.NphoMatched[rKey][i] = n4;

                f1[i] = (nLead>0)? n1/nLead : 0.0;
                f2[i] = (nLead>0)? n2/nLead : 0.0;
                f3[i] = (nLead>0)? n3/nLead : 0.0;
                f4[i] = (nLead>0)? n4/nLead : 0.0;

                totAll += nLead;
                tot1 += n1; tot2 += n2; tot3 += n3; tot4 += n4;
              }

                // Fractions plot (existing)
                {
                  TCanvas c("c_frac","c_frac",900,700);
                  ApplyCanvasMargins1D(c);

                  TGraph g1(kNPtBins, &x[0], &f1[0]);
                  TGraph g2(kNPtBins, &x[0], &f2[0]);
                  TGraph g3(kNPtBins, &x[0], &f3[0]);
                  TGraph g4(kNPtBins, &x[0], &f4[0]);

                  g1.SetLineWidth(2); g2.SetLineWidth(2); g3.SetLineWidth(2); g4.SetLineWidth(2);
                  g1.SetMarkerStyle(20); g2.SetMarkerStyle(21); g3.SetMarkerStyle(22); g4.SetMarkerStyle(24);
                  g1.SetLineColor(1); g2.SetLineColor(2); g3.SetLineColor(4); g4.SetLineColor(6);
                  g1.SetMarkerColor(1); g2.SetMarkerColor(2); g3.SetMarkerColor(4); g4.SetMarkerColor(6);

                  TMultiGraph mg;
                  mg.Add(&g1, "LP"); mg.Add(&g2, "LP"); mg.Add(&g3, "LP"); mg.Add(&g4, "LP");
                  mg.Draw("A");
                  mg.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                  mg.GetYaxis()->SetTitle("Fraction");
                  mg.SetMinimum(0.0); mg.SetMaximum(1.05);

                  TLegend leg(0.55,0.70,0.92,0.90);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.030);
                  leg.AddEntry(&g1, "NoJetPt", "lp");
                  leg.AddEntry(&g2, "NoJetEta", "lp");
                  leg.AddEntry(&g3, "NoBackToBack", "lp");
                  leg.AddEntry(&g4, "Matched", "lp");
                  leg.Draw();

                  DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                  DrawLatexLines(0.14,0.78, {string("Match status fractions vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);
                  SaveCanvas(c, JoinPath(D.dirProj, "match_status_fractions_vs_pTgamma.png"));
                }

                // 3x3 bar-table (per pT bin) of status fractions, saved into MatchQA/<rKey>/
                {
                  // Bin labels (Y categories in the original TH2):
                  // 1 = NoJetPt, 2 = NoJetEta, 3 = NoBackToBack, 4 = Matched
                  const char* xLabels[4] = {
                    "No jet passes p_{T}^{min}",
                    "Jet fails fiducial |#eta| < (1.1 - R)",
                    "Jets not back-to-back with #gamma",
                    "Valid recoil jet found"
                  };

                  // Solid colors per bar (match your request):
                  // 1 red, 2 blue, 3 dark green, 4 dark orange
                  const int binColors[4] = {
                    kRed+1,
                    kBlue+1,
                    kGreen+2,
                    kOrange+7
                  };

                  TCanvas ctbl("c_matchStatusBars","c_matchStatusBars",1500,1050);
                  ctbl.Divide(3,3,0.001,0.001);

                  // Keep all drawn hist objects alive until after SaveCanvas
                  std::vector<TObject*> keep;
                  keep.reserve(kNPtBins * 6); // axis + 4 bars (+ a little slack) per pad

                  for (int i = 0; i < kNPtBins; ++i)
                  {
                    const PtBin& b = PtBins()[i];
                    ctbl.cd(i+1);

                    gPad->SetTicks(1,1);
                    gPad->SetLeftMargin(0.16);
                    gPad->SetRightMargin(0.05);
                    gPad->SetTopMargin(0.14);
                    gPad->SetBottomMargin(0.40);

                    const double vals[4] = { f1[i], f2[i], f3[i], f4[i] };
                    const double yMaxPlot = 1.15; // give headroom for numeric labels

                    // Axis-only frame (no shading)
                    TH1F* hAxis = new TH1F(
                      TString::Format("h_matchStatusFracAxis_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data(),
                      "",
                      4, 0.5, 4.5
                    );
                    hAxis->SetDirectory(nullptr);
                    hAxis->SetStats(0);
                    hAxis->SetMinimum(0.0);
                    hAxis->SetMaximum(yMaxPlot);

                      for (int ib = 1; ib <= 4; ++ib)
                      {
                        hAxis->GetXaxis()->SetBinLabel(ib, TString::Format("%d", ib).Data());
                      }

                    hAxis->GetYaxis()->SetTitle("Fraction (within p_{T}^{#gamma} bin)");
                    hAxis->GetXaxis()->SetTitle("");
                    hAxis->GetXaxis()->LabelsOption("v");

                    hAxis->GetYaxis()->SetTitleSize(0.07);
                    hAxis->GetYaxis()->SetLabelSize(0.06);
                    hAxis->GetYaxis()->SetTitleOffset(1.10);

                    hAxis->GetXaxis()->SetLabelSize(0.09);
                    hAxis->GetXaxis()->SetLabelOffset(0.015);

                    hAxis->SetLineColor(kBlack);
                    hAxis->SetLineWidth(2);
                    hAxis->SetFillStyle(0);
                    hAxis->Draw("hist");

                    // Draw solid 2D bars (no BAR2 shading) — one per bin, colored
                    for (int ib = 1; ib <= 4; ++ib)
                    {
                      TH1F* hb = new TH1F(
                        TString::Format("h_matchStatusFracBar_%s_%s_%d_b%d",
                          ds.label.c_str(), rKey.c_str(), i, ib).Data(),
                        "",
                        4, 0.5, 4.5
                      );
                      hb->SetDirectory(nullptr);
                      hb->SetStats(0);

                      hb->SetBinContent(ib, vals[ib-1]);

                      hb->SetFillStyle(1001);
                      hb->SetFillColor(binColors[ib-1]);
                      hb->SetLineColor(kBlack);
                      hb->SetLineWidth(2);

                      hb->SetBarWidth(0.90);
                      hb->SetBarOffset(0.05);

                      hb->Draw("BAR SAME");

                      keep.push_back(hb);
                    }

                    // Print numeric value above each bar (like table3x3_preselectionFails.png)
                    TLatex t;
                    t.SetTextFont(42);
                    t.SetTextAlign(22);   // centered
                    t.SetTextSize(0.075); // tuned for small pads

                    for (int ib = 1; ib <= 4; ++ib)
                    {
                      const double y = vals[ib-1];
                      if (y <= 0.0) continue;

                      const double xC = hAxis->GetXaxis()->GetBinCenter(ib);
                      const double yText = std::min(y + 0.03*yMaxPlot, 0.95*yMaxPlot);
                      t.DrawLatex(xC, yText, TString::Format("%.2f", y).Data());
                    }

                    // pT-bin label at top of each pad
                    DrawLatexLines(0.16, 0.90, {b.label}, 0.085, 0.10);

                    keep.push_back(hAxis);
                  }

                  DrawLatexLines(0.12,0.98,
                    { "MatchQA: per-bin status fractions (bars)", rKey + RLabel(rKey) },
                    0.030, 0.040
                  );

                  // per-jet cutflow status fractions (NOT event-level),
                  // using h_jetcutflow_status_vs_pTgamma_<rKey>
                  {
                      TH2* hJet = GetObj<TH2>(ds, "h_jetcutflow_status_vs_pTgamma_" + rKey, false, false, false);
                      if (!hJet)
                      {
                        notes.push_back("Missing h_jetcutflow_status_vs_pTgamma_" + rKey);
                      }
                      else
                      {
                        vector<double> jf1(kNPtBins, 0.0), jf2(kNPtBins, 0.0), jf3(kNPtBins, 0.0), jf4(kNPtBins, 0.0);

                        for (int i = 0; i < kNPtBins; ++i)
                        {
                          const int xbin = i+1;
                          const double n1 = hJet->GetBinContent(xbin, 1);
                          const double n2 = hJet->GetBinContent(xbin, 2);
                          const double n3 = hJet->GetBinContent(xbin, 3);
                          const double n4 = hJet->GetBinContent(xbin, 4);
                          const double nTot = n1+n2+n3+n4;

                          jf1[i] = (nTot>0)? n1/nTot : 0.0;
                          jf2[i] = (nTot>0)? n2/nTot : 0.0;
                          jf3[i] = (nTot>0)? n3/nTot : 0.0;
                          jf4[i] = (nTot>0)? n4/nTot : 0.0;
                        }

                        // Reuse your same 3x3 colored bar-table style
                        const int binColors[4] = { kRed+1, kBlue+1, kGreen+2, kOrange+7 };

                        TCanvas ctblJ("c_jetCutflowBars","c_jetCutflowBars",1500,1050);
                        ctblJ.Divide(3,3,0.001,0.001);

                        std::vector<TObject*> keepJ;
                        keepJ.reserve(kNPtBins * 6);

                        for (int i = 0; i < kNPtBins; ++i)
                        {
                          const PtBin& b = PtBins()[i];
                          ctblJ.cd(i+1);

                          gPad->SetTicks(1,1);
                          gPad->SetLeftMargin(0.16);
                          gPad->SetRightMargin(0.05);
                          gPad->SetTopMargin(0.14);
                          gPad->SetBottomMargin(0.25);

                          const double vals[4] = { jf1[i], jf2[i], jf3[i], jf4[i] };
                          const double yMaxPlot = 1.15;

                          TH1F* hAxis = new TH1F(
                            TString::Format("h_jetCutflowAxis_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data(),
                            "",
                            4, 0.5, 4.5
                          );
                          hAxis->SetDirectory(nullptr);
                          hAxis->SetStats(0);
                          hAxis->SetMinimum(0.0);
                          hAxis->SetMaximum(yMaxPlot);

                          for (int ib = 1; ib <= 4; ++ib)
                            hAxis->GetXaxis()->SetBinLabel(ib, TString::Format("%d", ib).Data());

                          hAxis->GetYaxis()->SetTitle("Fraction of jets (within p_{T}^{#gamma} bin)");
                          hAxis->GetXaxis()->SetTitle("");

                          hAxis->GetYaxis()->SetTitleSize(0.07);
                          hAxis->GetYaxis()->SetLabelSize(0.06);
                          hAxis->GetYaxis()->SetTitleOffset(1.10);

                          hAxis->GetXaxis()->SetLabelSize(0.11);
                          hAxis->GetXaxis()->SetLabelOffset(0.010);

                          hAxis->SetLineColor(kBlack);
                          hAxis->SetLineWidth(2);
                          hAxis->SetFillStyle(0);
                          hAxis->Draw("hist");

                          for (int ib = 1; ib <= 4; ++ib)
                          {
                            TH1F* hb = new TH1F(
                              TString::Format("h_jetCutflowBar_%s_%s_%d_b%d",
                                ds.label.c_str(), rKey.c_str(), i, ib).Data(),
                              "",
                              4, 0.5, 4.5
                            );
                            hb->SetDirectory(nullptr);
                            hb->SetStats(0);

                            hb->SetBinContent(ib, vals[ib-1]);

                            hb->SetFillStyle(1001);
                            hb->SetFillColor(binColors[ib-1]);
                            hb->SetLineColor(kBlack);
                            hb->SetLineWidth(2);

                            hb->SetBarWidth(0.90);
                            hb->SetBarOffset(0.05);

                            hb->Draw("BAR SAME");
                            keepJ.push_back(hb);
                          }

                          TLatex t;
                          t.SetTextFont(42);
                          t.SetTextAlign(22);
                          t.SetTextSize(0.075);

                          for (int ib = 1; ib <= 4; ++ib)
                          {
                            const double y = vals[ib-1];
                            const double xC = hAxis->GetXaxis()->GetBinCenter(ib);
                            const double yText = std::min(y + 0.03*yMaxPlot, 0.95*yMaxPlot);
                            t.DrawLatex(xC, yText, TString::Format("%.2f", y).Data());
                          }

                          DrawLatexLines(0.16, 0.90, {b.label}, 0.085, 0.10);
                          keepJ.push_back(hAxis);
                        }

                        DrawLatexLines(0.12,0.98,
                          { "MatchQA: jet cutflow status fractions (per jet; not iso/tight-conditioned)",
                            rKey + RLabel(rKey),
                            "Bin map: 1=pT fail, 2=|#eta| fail, 3=|#Delta#phi| fail, 4=pass all" },
                          0.030, 0.040
                        );

                        SaveCanvas(ctblJ, JoinPath(D.rOut, "jetcutflow_status_fraction_bars_3x3.png"));

                        for (auto* obj : keepJ) delete obj;
                      }
                    }

                  for (auto* obj : keep) delete obj;
              }

              // Efficiency-style summaries (your added functionality)
              // Definitions:
              //   f_pT  = (N2+N3+N4)/Ntot
              //   f_fid = (N3+N4)/Ntot
              //   f_b2b = (N4)/Ntot
              //   P(fid|pt)  = (N3+N4)/(N2+N3+N4)
              //   P(b2b|fid) = (N4)/(N3+N4)
              vector<double> fPt(kNPtBins, 0.0), fFid(kNPtBins, 0.0), fB2B(kNPtBins, 0.0);
              vector<double> pFidGivenPt(kNPtBins, 0.0), pB2BGivenFid(kNPtBins, 0.0);

              for (int i = 0; i < kNPtBins; ++i)
              {
                const int xbin = i+1;
                const double n1 = hStatus->GetBinContent(xbin, 1);
                const double n2 = hStatus->GetBinContent(xbin, 2);
                const double n3 = hStatus->GetBinContent(xbin, 3);
                const double n4 = hStatus->GetBinContent(xbin, 4);

                const double nTot     = n1 + n2 + n3 + n4;
                const double nPassPt  = n2 + n3 + n4;
                const double nPassFid = n3 + n4;

                fPt[i]          = SafeDivide(nPassPt,  nTot,     0.0);
                fFid[i]         = SafeDivide(nPassFid, nTot,     0.0);
                fB2B[i]         = SafeDivide(n4,       nTot,     0.0);
                pFidGivenPt[i]  = SafeDivide(nPassFid, nPassPt,  0.0);
                pB2BGivenFid[i] = SafeDivide(n4,       nPassFid, 0.0);
              }

              const string dirEff = JoinPath(D.rOut, "efficiencies");
              EnsureDir(dirEff);

              // (A) Selection efficiencies vs pTgamma
              {
                TCanvas c2("c_effA","c_effA",900,700);
                ApplyCanvasMargins1D(c2);

                TGraph gPt(kNPtBins, &x[0], &fPt[0]);
                TGraph gFid(kNPtBins, &x[0], &fFid[0]);
                TGraph gB2B(kNPtBins, &x[0], &fB2B[0]);

                gPt.SetLineWidth(2);   gFid.SetLineWidth(2);   gB2B.SetLineWidth(2);
                gPt.SetMarkerStyle(20); gFid.SetMarkerStyle(21); gB2B.SetMarkerStyle(24);
                gPt.SetLineColor(1);   gFid.SetLineColor(2);   gB2B.SetLineColor(4);
                gPt.SetMarkerColor(1); gFid.SetMarkerColor(2); gB2B.SetMarkerColor(4);

                TMultiGraph mg2;
                mg2.Add(&gPt,  "LP");
                mg2.Add(&gFid, "LP");
                mg2.Add(&gB2B, "LP");

                mg2.Draw("A");
                mg2.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                mg2.GetYaxis()->SetTitle("Efficiency / fraction");
                mg2.SetMinimum(0.0);
                mg2.SetMaximum(1.05);

                TLegend leg2(0.55,0.72,0.92,0.90);
                leg2.SetTextFont(42);
                leg2.SetTextSize(0.030);
                leg2.AddEntry(&gPt,  "f_{pT}: jet exists above p_{T} cut", "lp");
                leg2.AddEntry(&gFid, "f_{fid}: fiducial jet exists",       "lp");
                leg2.AddEntry(&gB2B, "f_{b2b}: back-to-back recoil exists", "lp");
                leg2.Draw();

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, { "Photon-conditioned selection efficiencies", rKey }, 0.030, 0.040);

                SaveCanvas(c2, JoinPath(dirEff, "match_efficiencies_selection_vs_pTgamma.png"));
              }

              // (B) Conditional diagnostics vs pTgamma
              {
                TCanvas c3("c_effB","c_effB",900,700);
                ApplyCanvasMargins1D(c3);

                TGraph gA(kNPtBins, &x[0], &pFidGivenPt[0]);
                TGraph gB(kNPtBins, &x[0], &pB2BGivenFid[0]);

                gA.SetLineWidth(2); gB.SetLineWidth(2);
                gA.SetMarkerStyle(20); gB.SetMarkerStyle(24);
                gA.SetLineColor(1); gB.SetLineColor(2);
                gA.SetMarkerColor(1); gB.SetMarkerColor(2);

                TMultiGraph mg3;
                mg3.Add(&gA, "LP");
                mg3.Add(&gB, "LP");

                mg3.Draw("A");
                mg3.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                mg3.GetYaxis()->SetTitle("Conditional probability");
                mg3.SetMinimum(0.0);
                mg3.SetMaximum(1.05);

                TLegend leg3(0.55,0.75,0.92,0.90);
                leg3.SetTextFont(42);
                leg3.SetTextSize(0.030);
                leg3.AddEntry(&gA, "P(fid | pass p_{T})",  "lp");
                leg3.AddEntry(&gB, "P(b2b | fid)",         "lp");
                leg3.Draw();

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, { "Conditional matching diagnostics", rKey }, 0.030, 0.040);

                SaveCanvas(c3, JoinPath(dirEff, "match_efficiencies_conditionals_vs_pTgamma.png"));
              }

              // Terminal table + summary.txt (your added functionality)
              {
                cout << ANSI_BOLD_CYN
                     << "\n[MatchStatus table] " << ds.label << "  rKey=" << rKey
                     << " (R=" << std::fixed << std::setprecision(1) << RFromKey(rKey) << ")\n"
                     << ANSI_RESET;

                const int wPt   = 10;
                const int wN    = 12;
                const int wFrac = 10;

                cout << std::left << std::setw(wPt) << "pTbin"
                     << std::right
                     << std::setw(wN)    << "NoJetPt"
                     << std::setw(wN)    << "NoJetEta"
                     << std::setw(wN)    << "NoB2B"
                     << std::setw(wN)    << "Matched"
                     << std::setw(wN)    << "Total"
                     << "  |  "
                     << std::setw(wFrac) << "fMatch"
                     << std::setw(wFrac) << "fNoPt"
                     << std::setw(wFrac) << "fNoEta"
                     << std::setw(wFrac) << "fNoB2B"
                     << "  |  "
                     << std::setw(wFrac) << "P(fid|pt)"
                     << std::setw(wFrac) << "P(b2b|fid)"
                     << "\n";

                cout << string(wPt + 5*wN + 3 + 4*wFrac + 3 + 2*wFrac, '-') << "\n";

                // Per pT bin
                for (int i = 0; i < kNPtBins; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  const int xbin = i + 1;

                  const double n1 = hStatus->GetBinContent(xbin, 1);
                  const double n2 = hStatus->GetBinContent(xbin, 2);
                  const double n3 = hStatus->GetBinContent(xbin, 3);
                  const double n4 = hStatus->GetBinContent(xbin, 4);

                  const double nTot    = n1 + n2 + n3 + n4;
                  const double nPassPt = n2 + n3 + n4;
                  const double nPassFid = n3 + n4;

                  const double fMatch = SafeDivide(n4, nTot, 0.0);
                  const double fNoPt  = SafeDivide(n1, nTot, 0.0);
                  const double fNoEta = SafeDivide(n2, nTot, 0.0);
                  const double fNoB2B = SafeDivide(n3, nTot, 0.0);

                  const double p_fid_given_pt  = SafeDivide(nPassFid, nPassPt, 0.0);
                  const double p_b2b_given_fid = SafeDivide(n4, nPassFid, 0.0);

                  cout << std::left << std::setw(wPt) << b.label
                       << std::right
                       << std::setw(wN) << std::fixed << std::setprecision(0) << n1
                       << std::setw(wN) << n2
                       << std::setw(wN) << n3
                       << std::setw(wN) << n4
                       << std::setw(wN) << nTot
                       << "  |  "
                       << std::setw(wFrac) << std::fixed << std::setprecision(4) << fMatch
                       << std::setw(wFrac) << fNoPt
                       << std::setw(wFrac) << fNoEta
                       << std::setw(wFrac) << fNoB2B
                       << "  |  "
                       << std::setw(wFrac) << p_fid_given_pt
                       << std::setw(wFrac) << p_b2b_given_fid
                       << "\n";
                }

                // Integrated summary
                const double totAll = (mc.NphoLead[rKey][0] + mc.NphoLead[rKey][1] + mc.NphoLead[rKey][2] +
                                       mc.NphoLead[rKey][3] + mc.NphoLead[rKey][4] + mc.NphoLead[rKey][5] +
                                       mc.NphoLead[rKey][6] + mc.NphoLead[rKey][7] + mc.NphoLead[rKey][8]);

                // Recompute integrated totals directly from TH2 (exact)
                double I1=0, I2=0, I3=0, I4=0, IAll=0;
                for (int i = 0; i < kNPtBins; ++i)
                {
                  const int xbin = i+1;
                  const double n1 = hStatus->GetBinContent(xbin, 1);
                  const double n2 = hStatus->GetBinContent(xbin, 2);
                  const double n3 = hStatus->GetBinContent(xbin, 3);
                  const double n4 = hStatus->GetBinContent(xbin, 4);
                  const double nTot = n1+n2+n3+n4;
                  I1 += n1; I2 += n2; I3 += n3; I4 += n4; IAll += nTot;
                }

                const double fMatchAll = SafeDivide(I4, IAll, 0.0);
                const double fNoPtAll  = SafeDivide(I1, IAll, 0.0);
                const double fNoEtaAll = SafeDivide(I2, IAll, 0.0);
                const double fNoB2BAll = SafeDivide(I3, IAll, 0.0);

                const double passPtAll  = I2 + I3 + I4;
                const double passFidAll = I3 + I4;

                const double p_fid_given_pt_all  = SafeDivide(passFidAll, passPtAll, 0.0);
                const double p_b2b_given_fid_all = SafeDivide(I4, passFidAll, 0.0);

                cout << string(10 + 5*12 + 3 + 4*10 + 3 + 2*10, '-') << "\n";
                cout << ANSI_BOLD_YEL
                     << "Integrated: N=" << std::fixed << std::setprecision(0) << IAll
                     << "  fMatch=" << std::setprecision(4) << fMatchAll
                     << "  fNoPt=" << fNoPtAll
                     << "  fNoEta=" << fNoEtaAll
                     << "  fNoB2B=" << fNoB2BAll
                     << "  P(fid|pt)=" << p_fid_given_pt_all
                     << "  P(b2b|fid)=" << p_b2b_given_fid_all
                     << ANSI_RESET << "\n";

                vector<string> s;
                s.push_back(string("MatchQA summary (") + ds.label + ")");
                s.push_back(string("rKey: ") + rKey + TString::Format("  R=%.1f", RFromKey(rKey)).Data());
                s.push_back("");
                s.push_back(TString::Format("Total leading photons: %.0f", IAll).Data());
                s.push_back("Integrated fractions:");
                s.push_back(TString::Format("  f_NoJetPt      = %.6f", SafeDivide(I1, IAll, 0.0)).Data());
                s.push_back(TString::Format("  f_NoJetEta     = %.6f", SafeDivide(I2, IAll, 0.0)).Data());
                s.push_back(TString::Format("  f_NoBackToBack = %.6f", SafeDivide(I3, IAll, 0.0)).Data());
                s.push_back(TString::Format("  f_Matched      = %.6f", SafeDivide(I4, IAll, 0.0)).Data());
                s.push_back("");
                s.push_back("Conditional diagnostics (integrated):");
                s.push_back(TString::Format("  P(fiducial | pass pT)  = %.6f", p_fid_given_pt_all).Data());
                s.push_back(TString::Format("  P(back-to-back | fid)  = %.6f", p_b2b_given_fid_all).Data());
                s.push_back("");
                s.push_back("Per pT bin:");
                for (int i = 0; i < kNPtBins; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  const int xbin = i + 1;

                  const double n1 = hStatus->GetBinContent(xbin, 1);
                  const double n2 = hStatus->GetBinContent(xbin, 2);
                  const double n3 = hStatus->GetBinContent(xbin, 3);
                  const double n4 = hStatus->GetBinContent(xbin, 4);
                  const double nTot = n1 + n2 + n3 + n4;

                  const double nPassPt = n2 + n3 + n4;
                  const double nPassFid = n3 + n4;

                  const double p_fid_given_pt  = SafeDivide(nPassFid, nPassPt, 0.0);
                  const double p_b2b_given_fid = SafeDivide(n4, nPassFid, 0.0);

                  s.push_back(TString::Format(
                    "  %s  N=%.0f  NoPt=%.0f  NoEta=%.0f  NoB2B=%.0f  Match=%.0f  fMatch=%.6f  P(fid|pt)=%.6f  P(b2b|fid)=%.6f",
                    b.label.c_str(), nTot, n1, n2, n3, n4,
                    SafeDivide(n4, nTot, 0.0),
                    p_fid_given_pt,
                    p_b2b_given_fid
                  ).Data());
                }

                if (!notes.empty())
                {
                  s.push_back("");
                  s.push_back("NOTES:");
                  for (const auto& n : notes) s.push_back(string("  - ") + n);
                }

                WriteTextFile(JoinPath(D.rOut, "summary.txt"), s);
              }

              return hStatus;
            };

            // -------------------------------------------------------------------------
            // Helper: maxdphi QA (existing behavior)
            // -------------------------------------------------------------------------
            auto HandleMaxDphi =
              [&](const string& rKey, const MatchDirs& D)
            {
              if (TH2* hMax = GetObj<TH2>(ds, "h_match_maxdphi_vs_pTgamma_" + rKey, true, true, true))
              {
                TH2* hc = CloneTH2(hMax, "maxdphi_clone");
                DrawAndSaveTH2_Common(ds, hc,
                  JoinPath(D.dir2D, "match_maxdphi_vs_pTgamma.png"),
                  "p_{T}^{#gamma} [GeV]", "max|#Delta#phi| [rad]", "Counts",
                  {string("max|#Delta#phi(#gamma,jet)| vs p_{T}^{#gamma}"), rKey}, false);
                delete hc;

                vector<double> x(kNPtBins, 0.0), y(kNPtBins, 0.0);
                for (int i = 0; i < kNPtBins; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  x[i] = 0.5*(b.lo + b.hi);
                  TH1D* proj = hMax->ProjectionY("tmp_maxdphi_py", i+1, i+1);
                  y[i] = proj ? proj->GetMean() : 0.0;
                  if (proj) delete proj;
                }

                TCanvas c("c_meanmax","c_meanmax",900,700);
                ApplyCanvasMargins1D(c);
                TGraph g(kNPtBins, &x[0], &y[0]);
                g.SetLineWidth(2);
                g.SetMarkerStyle(20);
                g.Draw("ALP");
                g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g.GetYaxis()->SetTitle("<max|#Delta#phi|> [rad]");
                g.SetMinimum(0.0);
                g.SetMaximum(TMath::Pi());

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {string("Mean max|#Delta#phi| vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);

                SaveCanvas(c, JoinPath(D.dirProj, "mean_maxdphi_vs_pTgamma.png"));
              }
            };

            // -------------------------------------------------------------------------
            // Helper: profile nRecoil jets vs pTgamma (existing behavior)
            // -------------------------------------------------------------------------
            auto HandleNRecoilProfile =
              [&](const string& rKey, const MatchDirs& D)
            {
              if (TProfile* pN = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_" + rKey, true, true, true))
              {
                TProfile* pc = (TProfile*)pN->Clone("p_clone");
                pc->SetDirectory(nullptr);
                DrawAndSaveTH1_Common(ds, (TH1*)pc,
                  JoinPath(D.dirProj, "profile_nRecoilJets_vs_pTgamma.png"),
                  "p_{T}^{#gamma} [GeV]", "<N_{recoil jets}>",
                  {string("Mean recoil-jet multiplicity vs p_{T}^{#gamma}"), rKey}, false);
                delete pc;
              }
            };

            // -------------------------------------------------------------------------
            // Helper: additional Δφ matching QA suite (your large added block, organized)
            // -------------------------------------------------------------------------
            auto HandleDeltaPhiSuite =
              [&](const string& rKey, const MatchDirs& D)
            {
              // Uses:
              //  - h_match_dphi_vs_pTgamma_rKey    : |Δphi(γ, recoilJet1)| for matched events only
              //  - h_match_maxdphi_vs_pTgamma_rKey : max|Δphi(γ, fid jet)| for all events (no-fid-jet -> underflow)
              TH2* hDphi = GetObj<TH2>(ds, "h_match_dphi_vs_pTgamma_" + rKey, true, true, true);
              if (!hDphi) return;

              // Keep existing 2D plot (matched recoilJet1 only)
              {
                TH2* hc = CloneTH2(hDphi,
                  TString::Format("dphi_clone_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                );

                DrawAndSaveTH2_Common(ds, hc,
                  JoinPath(D.dir2D, "match_dphi_vs_pTgamma.png"),
                  "p_{T}^{#gamma} [GeV]", "|#Delta#phi| [rad]", "Counts",
                  {string("|#Delta#phi(#gamma,recoilJet1)| vs p_{T}^{#gamma} (matched only)"), rKey},
                  false);

                delete hc;
              }

              // Optional companion TH2: max|Δphi| over fid jets (filled for all events)
              TH2* hMax = GetObj<TH2>(ds, "h_match_maxdphi_vs_pTgamma_" + rKey, false, false, false);

              // Output directories (clean + organized)
              const string dirDphi = JoinPath(D.dirProj, "deltaPhi");
              const string dirJet1 = JoinPath(dirDphi, "recoilJet1");
              const string dirMax  = JoinPath(dirDphi, "maxDphi");
              const string dirOv   = JoinPath(dirDphi, "overlays");

              EnsureDir(dirDphi);
              EnsureDir(dirJet1);
              EnsureDir(dirMax);
              EnsureDir(dirOv);

              for (const auto& b : PtBins())
              {
                EnsureDir(JoinPath(dirJet1, b.folder));
                EnsureDir(JoinPath(dirMax,  b.folder));
                EnsureDir(JoinPath(dirOv,   b.folder));
              }

              // Helpers: visible-area normalization + safe projections with unique names
              auto NormalizeVisible = [](TH1* h)
              {
                if (!h) return;
                const int nb = h->GetNbinsX();
                const double integral = h->Integral(1, nb); // exclude under/overflow
                if (integral > 0.0) h->Scale(1.0 / integral);
              };

              auto ProjY_AtXbin = [&](TH2* h2, int xbin, const string& newName)->TH1D*
              {
                if (!h2) return nullptr;
                TH1D* p = h2->ProjectionY(newName.c_str(), xbin, xbin);
                if (p) p->SetDirectory(nullptr);
                return p;
              };

              auto ProjY_AllX = [&](TH2* h2, const string& newName)->TH1D*
              {
                if (!h2) return nullptr;
                const int nx = h2->GetXaxis()->GetNbins();
                TH1D* p = h2->ProjectionY(newName.c_str(), 1, nx);
                if (p) p->SetDirectory(nullptr);
                return p;
              };

              const int nPtAxis = hDphi->GetXaxis()->GetNbins();
              const int nPtUse  = std::min(nPtAxis, kNPtBins);

              if (nPtAxis != kNPtBins)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] " << ds.label << " " << rKey
                     << ": h_match_dphi_vs_pTgamma has nPtBins=" << nPtAxis
                     << " but offline expects kNPtBins=" << kNPtBins
                     << " (will plot first " << nPtUse << " bins only)\n"
                     << ANSI_RESET;
              }

              // pT-INDEPENDENT: recoilJet1 |Δphi| integrated
              {
                TH1D* pAll = ProjY_AllX(
                  hDphi,
                  TString::Format("p_dphiJet1_all_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                );

                if (pAll)
                {
                  vector<string> lines = {
                    "Matched recoilJet1: |#Delta#phi(#gamma,jet1)|",
                    TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), RFromKey(rKey)).Data(),
                    "Integrated over all p_{T}^{#gamma} bins",
                    "Filled only when recoilJet1 exists (matchStatus=Matched)"
                  };

                  DrawAndSaveTH1_Common(ds, pAll,
                    JoinPath(dirJet1, "dphi_jet1_integrated_counts.png"),
                    "|#Delta#phi| [rad]", "Counts", lines, false);

                  TH1* pShape = CloneTH1(pAll,
                    TString::Format("p_dphiJet1_allShape_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                  );
                  NormalizeVisible(pShape);

                  DrawAndSaveTH1_Common(ds, pShape,
                    JoinPath(dirJet1, "dphi_jet1_integrated_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", lines, false);

                  delete pShape;
                  delete pAll;
                }
              }

              // pT-INDEPENDENT: max|Δphi| integrated
              TH1* maxShapeForIntegratedOverlay = nullptr;

              if (hMax)
              {
                TH1D* pAll = ProjY_AllX(
                  hMax,
                  TString::Format("p_maxDphi_all_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                );

                if (pAll)
                {
                  vector<string> lines = {
                    "max|#Delta#phi(#gamma,jet)| over fid jets (pass jet p_{T}+|#eta|)",
                    TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), RFromKey(rKey)).Data(),
                    "Integrated over all p_{T}^{#gamma} bins",
                    "NOTE: events with no fid jet fill underflow (not visible)"
                  };

                  DrawAndSaveTH1_Common(ds, pAll,
                    JoinPath(dirMax, "maxDphi_integrated_counts.png"),
                    "max|#Delta#phi| [rad]", "Counts", lines, false);

                  TH1* pShape = CloneTH1(pAll,
                    TString::Format("p_maxDphi_allShape_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                  );
                  NormalizeVisible(pShape);

                  DrawAndSaveTH1_Common(ds, pShape,
                    JoinPath(dirMax, "maxDphi_integrated_shape.png"),
                    "max|#Delta#phi| [rad]", "A.U.", lines, false);

                  maxShapeForIntegratedOverlay = CloneTH1(pShape,
                    TString::Format("p_maxDphi_allShape_forOv_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                  );

                  delete pShape;
                  delete pAll;
                }
              }

              // pT-INDEPENDENT overlay (shape): recoilJet1 vs max|Δphi|
              if (hMax && maxShapeForIntegratedOverlay)
              {
                TH1D* pJet1All = ProjY_AllX(
                  hDphi,
                  TString::Format("p_dphiJet1_all_forOv_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                );

                if (pJet1All)
                {
                  TH1* jet1Shape = CloneTH1(pJet1All,
                    TString::Format("p_dphiJet1_allShape_forOv_%s_%s", ds.label.c_str(), rKey.c_str()).Data()
                  );
                  NormalizeVisible(jet1Shape);

                  vector<string> extra = {
                    "Overlay (shape): recoilJet1 vs max|#Delta#phi|",
                    "If curves differ: leading recoil jet is not always the most back-to-back jet"
                  };

                  DrawOverlayTwoTH1(ds, jet1Shape, maxShapeForIntegratedOverlay,
                    "recoilJet1 (matched)", "max|#Delta#phi| (fid jets)",
                    JoinPath(dirOv, "overlay_dphiJet1_vs_maxDphi_integrated_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", extra, false);

                  delete jet1Shape;
                  delete pJet1All;
                }

                delete maxShapeForIntegratedOverlay;
                maxShapeForIntegratedOverlay = nullptr;
              }

              // PER-pT bin: projections + overlays
              for (int i = 0; i < nPtUse; ++i)
              {
                const PtBin& b = PtBins()[i];
                const int xbin = i + 1;

                TH1D* p1 = ProjY_AtXbin(
                  hDphi, xbin,
                  TString::Format("p_dphiJet1_%s_%s_%s", ds.label.c_str(), rKey.c_str(), b.folder.c_str()).Data()
                );

                if (p1)
                {
                  vector<string> lines = {
                    "Matched recoilJet1: |#Delta#phi(#gamma,jet1)|",
                    TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), RFromKey(rKey)).Data(),
                    TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data(),
                    "Filled only when recoilJet1 exists (matchStatus=Matched)"
                  };

                  DrawAndSaveTH1_Common(ds, p1,
                    JoinPath(dirJet1, b.folder + "/dphi_jet1_counts.png"),
                    "|#Delta#phi| [rad]", "Counts", lines, false);

                  TH1* p1Shape = CloneTH1(p1,
                    TString::Format("p_dphiJet1_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                  );
                  NormalizeVisible(p1Shape);

                  DrawAndSaveTH1_Common(ds, p1Shape,
                    JoinPath(dirJet1, b.folder + "/dphi_jet1_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", lines, false);

                  delete p1Shape;
                }

                TH1D* pM = nullptr;
                if (hMax)
                {
                  pM = ProjY_AtXbin(
                    hMax, xbin,
                    TString::Format("p_maxDphi_%s_%s_%s", ds.label.c_str(), rKey.c_str(), b.folder.c_str()).Data()
                  );

                  if (pM)
                  {
                    vector<string> lines = {
                      "max|#Delta#phi(#gamma,jet)| over fid jets (pass jet p_{T}+|#eta|)",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), RFromKey(rKey)).Data(),
                      TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data(),
                      "NOTE: no-fid-jet events fill underflow (not visible)"
                    };

                    DrawAndSaveTH1_Common(ds, pM,
                      JoinPath(dirMax, b.folder + "/maxDphi_counts.png"),
                      "max|#Delta#phi| [rad]", "Counts", lines, false);

                    TH1* pMShape = CloneTH1(pM,
                      TString::Format("p_maxDphi_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                    );
                    NormalizeVisible(pMShape);

                    DrawAndSaveTH1_Common(ds, pMShape,
                      JoinPath(dirMax, b.folder + "/maxDphi_shape.png"),
                      "max|#Delta#phi| [rad]", "A.U.", lines, false);

                    delete pMShape;
                  }
                }

                // Overlay (shape): recoilJet1 vs max|Δphi|
                if (p1 && pM)
                {
                  TH1* a = CloneTH1(p1,
                    TString::Format("ov_dphiJet1_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                  );
                  TH1* b2 = CloneTH1(pM,
                    TString::Format("ov_maxDphi_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                  );

                  NormalizeVisible(a);
                  NormalizeVisible(b2);

                  vector<string> extra = {
                    TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data(),
                    "Overlay (shape): recoilJet1 vs max|#Delta#phi|"
                  };

                  DrawOverlayTwoTH1(ds, a, b2,
                    "recoilJet1 (matched)", "max|#Delta#phi| (fid jets)",
                    JoinPath(dirOv, b.folder + "/overlay_dphiJet1_vs_maxDphi_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", extra, false);

                  delete a;
                  delete b2;
                }

                if (p1) delete p1;
                if (pM) delete pM;
              }

              // 3x3 tables (shape): jet1 dphi and maxDphi and overlay table
              auto Make3x3Table_ProjYShape =
                [&](TH2* h2,
                    const string& canvasTag,
                    const string& outPng,
                    const string& xTitle,
                    const vector<string>& headerLines,
                    bool logy)
              {
                if (!h2) return;

                TCanvas c(
                  TString::Format("c_tbl_dphi_%s_%s_%s", ds.label.c_str(), rKey.c_str(), canvasTag.c_str()).Data(),
                  "c_tbl_dphi", 1500, 900
                );
                c.Divide(3,2, 0.001, 0.001);

                vector<TH1*> keep;
                keep.reserve(kNPtBins);

                const int nPads = std::min(kNPtBins, nPtUse);
                for (int i = 0; i < nPads; ++i)
                {
                  c.cd(i+1);
                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.14);
                  gPad->SetTopMargin(0.10);
                  gPad->SetLogy(logy);

                  const PtBin& b = PtBins()[i];
                  const int xbin = i + 1;

                  TH1D* p = ProjY_AtXbin(
                    h2, xbin,
                    TString::Format("tbl_proj_%s_%s_%d_%s", ds.label.c_str(), rKey.c_str(), i, canvasTag.c_str()).Data()
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

                  NormalizeVisible(p);
                  p->SetLineWidth(2);
                  p->SetTitle("");
                  p->GetXaxis()->SetTitle(xTitle.c_str());
                  p->GetYaxis()->SetTitle("A.U.");
                  p->Draw("hist");

                  vector<string> lines = headerLines;
                  lines.push_back(TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                  DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);

                  keep.push_back(p);
                }

                SaveCanvas(c, outPng);

                for (auto* h : keep) delete h;
              };

              // jet1 dphi table
              {
                vector<string> hdr = {
                  "Matched recoilJet1: |#Delta#phi(#gamma,jet1)| (shape)",
                  TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), RFromKey(rKey)).Data()
                };
                Make3x3Table_ProjYShape(
                  hDphi,
                  "jet1",
                  JoinPath(dirJet1, "table3x3_dphi_jet1_shape.png"),
                  "|#Delta#phi| [rad]",
                  hdr,
                  false
                );
              }

              // maxDphi table
              if (hMax)
              {
                vector<string> hdr = {
                  "max|#Delta#phi(#gamma,jet)| over fid jets (shape)",
                  TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), RFromKey(rKey)).Data()
                };
                Make3x3Table_ProjYShape(
                  hMax,
                  "max",
                  JoinPath(dirMax, "table3x3_maxDphi_shape.png"),
                  "max|#Delta#phi| [rad]",
                  hdr,
                  false
                );
              }

              // overlay table jet1 vs max
              if (hMax)
              {
                TCanvas c(
                  TString::Format("c_tbl_dphiOv_%s_%s", ds.label.c_str(), rKey.c_str()).Data(),
                  "c_tbl_dphiOv", 1500, 1200
                );
                c.Divide(3,3, 0.001, 0.001);

                vector<TObject*> keep;
                keep.reserve(3 * kNPtBins);

                for (int i = 0; i < kNPtBins; ++i)
                {
                  c.cd(i+1);
                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.14);
                  gPad->SetTopMargin(0.10);

                  if (i >= nPtUse)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.20, 0.55, "EMPTY");
                    continue;
                  }

                  const PtBin& b = PtBins()[i];
                  const int xbin = i + 1;

                  TH1D* p1 = ProjY_AtXbin(
                    hDphi, xbin,
                    TString::Format("tblOv_p1_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                  );
                  TH1D* p2 = ProjY_AtXbin(
                    hMax, xbin,
                    TString::Format("tblOv_p2_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                  );

                  if (!p1 || !p2 || p1->GetEntries() <= 0.0 || p2->GetEntries() <= 0.0)
                  {
                    if (p1) delete p1;
                    if (p2) delete p2;
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.15, 0.55, "MISSING");
                    continue;
                  }

                  NormalizeVisible(p1);
                  NormalizeVisible(p2);

                  p1->SetLineWidth(2);
                  p2->SetLineWidth(2);
                  p1->SetLineColor(1);
                  p2->SetLineColor(2);

                  const double ymax = std::max(p1->GetMaximum(), p2->GetMaximum());
                  p1->SetMaximum(ymax * 1.25);

                  p1->SetTitle("");
                  p1->GetXaxis()->SetTitle("|#Delta#phi| [rad]");
                  p1->GetYaxis()->SetTitle("A.U.");

                  p1->Draw("hist");
                  p2->Draw("hist same");

                  TLegend* leg = new TLegend(0.55, 0.72, 0.92, 0.90);
                  leg->SetTextFont(42);
                  leg->SetTextSize(0.030);
                  leg->AddEntry(p1, "recoilJet1 (matched)", "l");
                  leg->AddEntry(p2, "max|#Delta#phi| (fid jets)", "l");
                  leg->Draw();

                  DrawLatexLines(0.16, 0.90,
                    {
                      TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data(),
                      "Overlay (shape)"
                    },
                    0.040, 0.050
                  );

                  keep.push_back(p1);
                  keep.push_back(p2);
                  keep.push_back(leg);
                }

                SaveCanvas(c, JoinPath(dirOv, "table3x3_overlay_dphiJet1_vs_maxDphi_shape.png"));

                for (auto* obj : keep) delete obj;
              }
            };

            // -------------------------------------------------------------------------
            // Helper: recoilIsLeading vs pTgamma + fraction plot
            // -------------------------------------------------------------------------
            auto HandleRecoilIsLeading =
              [&](const string& rKey, const MatchDirs& D)
            {
              if (TH2* hRL = GetObj<TH2>(ds, "h_recoilIsLeading_vs_pTgamma_" + rKey, true, true, true))
              {
                TH2* hc = CloneTH2(hRL, "recoilLead_clone");
                DrawAndSaveTH2_Common(ds, hc,
                  JoinPath(D.dir2D, "recoilIsLeading_vs_pTgamma.png"),
                  "p_{T}^{#gamma} [GeV]", "Status", "Counts",
                  {string("Recoil jet is leading? vs p_{T}^{#gamma}"), rKey}, false);
                delete hc;

                vector<double> x(kNPtBins, 0.0), y(kNPtBins, 0.0);
                for (int i = 0; i < kNPtBins; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  x[i] = 0.5*(b.lo + b.hi);
                  const double nNot = hRL->GetBinContent(i+1, 1);
                  const double nYes = hRL->GetBinContent(i+1, 2);
                  const double ntot = nNot + nYes;
                  y[i] = (ntot > 0.0) ? (nYes/ntot) : 0.0;
                }

                TCanvas c("c_flead","c_flead",900,700);
                ApplyCanvasMargins1D(c);
                TGraph g(kNPtBins, &x[0], &y[0]);
                g.SetLineWidth(2);
                g.SetMarkerStyle(20);
                g.Draw("ALP");
                g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g.GetYaxis()->SetTitle("f(recoilJet1 is leading)");
                g.SetMinimum(0.0); g.SetMaximum(1.05);

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {string("Fraction recoilJet1 is leading vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);

                SaveCanvas(c, JoinPath(D.dirProj, "fraction_recoilIsLeading_vs_pTgamma.png"));
              }
            };

            // -------------------------------------------------------------------------
            // MAIN LOOP: per rKey
            // -------------------------------------------------------------------------
            for (const auto& rKey : kRKeys)
            {
              const MatchDirs D = MakeDirsForRKey(rKey);

              // Keep notes for summary.txt
              vector<string> notes;

              // Match status (plots + cache fill + efficiencies + table + summary.txt)
              TH2* hStatus = HandleMatchStatus(rKey, D, notes);

              // Other QA blocks (independent)
              HandleMaxDphi(rKey, D);
              HandleNRecoilProfile(rKey, D);

              // Your extended Δφ suite (only runs if h_match_dphi exists)
              HandleDeltaPhiSuite(rKey, D);

              // recoil-is-leading
              HandleRecoilIsLeading(rKey, D);

              // If status histogram missing, we at least print a warning + optional note (keeps behavior safe)
              if (!hStatus && !notes.empty())
              {
                cout << ANSI_BOLD_YEL << "[WARN] MatchQA: missing status hist for " << ds.label << " " << rKey << ANSI_RESET << "\n";
              }
            }

            // -------------------------------------------------------------------------
            // OVERLAYS ACROSS rKeys (uses MatchCache filled from status hist)
            // -------------------------------------------------------------------------
            {
              const string overDir = JoinPath(baseOut, "overlays");
              EnsureDir(overDir);

                // x = pT bin centers; y = matched fraction; ey = binomial stat unc on fraction
                vector<double> x(kNPtBins, 0.0), ex(kNPtBins, 0.0);
                vector<double> y02(kNPtBins, 0.0), y04(kNPtBins, 0.0);
                vector<double> ey02(kNPtBins, 0.0), ey04(kNPtBins, 0.0);

                for (int i = 0; i < kNPtBins; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  x[i]  = 0.5*(b.lo + b.hi);
                  ex[i] = 0.0; // (optional) could use 0.5*(b.hi - b.lo) if you want x-bin-width bars

                  const double l02 = mc.NphoLead["r02"][i];
                  const double l04 = mc.NphoLead["r04"][i];
                  const double m02 = mc.NphoMatched["r02"][i];
                  const double m04 = mc.NphoMatched["r04"][i];

                  y02[i] = (l02 > 0.0) ? (m02/l02) : 0.0;
                  y04[i] = (l04 > 0.0) ? (m04/l04) : 0.0;

                  // Binomial proportion uncertainty: sqrt(p(1-p)/N)
                  ey02[i] = (l02 > 0.0) ? std::sqrt( y02[i]*(1.0 - y02[i]) / l02 ) : 0.0;
                  ey04[i] = (l04 > 0.0) ? std::sqrt( y04[i]*(1.0 - y04[i]) / l04 ) : 0.0;
                }

                // Matched fraction overlay (points only + stat error bars; no connecting lines)
                {
                  TCanvas c("c_ov_match","c_ov_match",900,700);
                  ApplyCanvasMargins1D(c);

                  TGraphErrors ge02(kNPtBins, &x[0], &y02[0], &ex[0], &ey02[0]);
                  TGraphErrors ge04(kNPtBins, &x[0], &y04[0], &ex[0], &ey04[0]);

                  // This controls the title that showed as "Graph" before
                  ge02.SetTitle("Recoil-Jet Match Probability (Reco, Radius-Dependent)");

                  ge02.SetLineWidth(2); ge04.SetLineWidth(2);
                  ge02.SetMarkerStyle(20); ge04.SetMarkerStyle(24);
                  ge02.SetLineColor(1);    ge04.SetLineColor(2);
                  ge02.SetMarkerColor(1);  ge04.SetMarkerColor(2);

                  // IMPORTANT: no "L" in the draw option -> no connecting lines
                  ge02.Draw("APE1");
                  ge02.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                  ge02.GetYaxis()->SetTitle("P(#geq1 jet passes as recoil for an accepeted reco #gamma)");
                  ge02.SetMinimum(0.0); ge02.SetMaximum(1.05);

                  ge04.Draw("PE1 same");

                  TLegend leg(0.62,0.22,0.92,0.45);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.033);
                  leg.AddEntry(&ge02, "R=0.2", "pe");
                  leg.AddEntry(&ge04, "R=0.4", "pe");
                  leg.Draw();

                  DrawLatexLines(0.14,0.85, DefaultHeaderLines(ds), 0.038, 0.045);

                  SaveCanvas(c, JoinPath(overDir, "overlay_matched_fraction_r02_vs_r04.png"));
                }


              // mean nRecoil jets overlay if profiles exist (existing)
              {
                TProfile* p02 = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_r02", true, true, true);
                TProfile* p04 = GetObj<TProfile>(ds, "p_nRecoilJets_vs_pTgamma_r04", true, true, true);
                if (p02 && p04)
                {
                  TProfile* a = (TProfile*)p02->Clone("p02c"); a->SetDirectory(nullptr);
                  TProfile* b = (TProfile*)p04->Clone("p04c"); b->SetDirectory(nullptr);

                  a->SetLineWidth(2); b->SetLineWidth(2);
                  a->SetMarkerStyle(20); b->SetMarkerStyle(24);
                  a->SetLineColor(1); b->SetLineColor(2);
                  a->SetMarkerColor(1); b->SetMarkerColor(2);

                  TCanvas c2("c_ov_nrecoil","c_ov_nrecoil",900,700);
                  ApplyCanvasMargins1D(c2);

                  a->SetTitle("");
                  a->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  a->GetYaxis()->SetTitle("<N_{recoil jets}>");
                  a->Draw("E1");
                  b->Draw("E1 same");

                  TLegend leg2(0.62,0.25,0.92,0.4);
                  leg2.SetTextFont(42);
                  leg2.SetTextSize(0.035);
                  leg2.AddEntry(a, "r02 (R=0.2)", "lp");
                  leg2.AddEntry(b, "r04 (R=0.4)", "lp");
                  leg2.Draw();

                  DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                  DrawLatexLines(0.14,0.78, {"Overlay: <N_{recoil jets}> vs p_{T}^{#gamma}"}, 0.030, 0.040);

                  SaveCanvas(c2, JoinPath(overDir, "overlay_mean_nRecoilJets_r02_vs_r04.png"));

                  delete a;
                  delete b;
                }
              }
           }
        }

        // =============================================================================
        // Section 5D: SelectedJetQA (jet1/jet2 distributions per pTgamma bin)
        // Uses:
        //   h_jet1Pt_rXX_pT_lo_hi
        //   h_jet1Eta_sel_rXX_pT_lo_hi
        //   h_jet1Phi_sel_rXX_pT_lo_hi
        //   h_jet1Mass_sel_rXX_pT_lo_hi
        //   and jet2 analogs (Pt/Phi/Eta_sel)
        // =============================================================================
        void RunSelectedJetQA(Dataset& ds)
        {
          cout << ANSI_BOLD_CYN << "\n==============================\n"
               << "[SECTION 5D] SelectedJetQA (jet1/jet2, per pT^{#gamma} bin) (" << ds.label << ")\n"
               << "==============================" << ANSI_RESET << "\n";

          string outDir = ds.isSim
              ? JoinPath(ds.outBase, "RecoilJetQA/SelectedJetQA")
              : JoinPath(ds.outBase, "baselineData/RecoilJetQA/SelectedJetQA");

          EnsureDir(outDir);
          for (const auto& rKey : kRKeys) EnsureDir(JoinPath(outDir, rKey));

          vector<string> common;
          common.push_back("Selected jets after #gamma-jet matching");
          common.push_back("Jets: p_{T}^{jet} #geq 10 GeV");

          for (const auto& rKey : kRKeys)
          {
            const double R = RFromKey(rKey);
            const double etaFidAbs = FidEtaAbsFromKey(rKey);

            const string rOut = JoinPath(outDir, rKey);
            EnsureDir(rOut);

            const string dirJet1 = JoinPath(rOut, "jet1");
            const string dirJet2 = JoinPath(rOut, "jet2");
            EnsureDir(dirJet1);
            EnsureDir(dirJet2);

            // Terminal summary header
            cout << ANSI_BOLD_CYN
                 << "\n[SelectedJetQA summary] " << ds.label << "  rKey=" << rKey
                 << " (R=" << std::fixed << std::setprecision(1) << R << ")\n"
                 << ANSI_RESET;

            const int wPt = 10;
            const int wN  = 12;
            const int wM  = 12;

            cout << std::left << std::setw(wPt) << "pTbin"
                 << std::right
                 << std::setw(wN) << "N(jet1)"
                 << std::setw(wN) << "N(jet2)"
                 << std::setw(wM) << "<pT1>"
                 << std::setw(wM) << "<pT2>"
                 << "\n";
            cout << string(wPt + 2*wN + 2*wM, '-') << "\n";

            vector<string> summaryLines;
            summaryLines.push_back(string("SelectedJetQA summary (") + ds.label + ")");
            summaryLines.push_back(string("rKey: ") + rKey + TString::Format("  R=%.1f", R).Data());
            summaryLines.push_back("");

            // Plot helpers
            auto plotVar = [&](const string& histBase, const string& subDir,
                               const string& outStem, const string& xTitle,
                               bool logy, bool shapeTable, bool drawFidLines)
            {
              const string sub = JoinPath(subDir, outStem);
              EnsureDir(sub);
              for (const auto& b : PtBins()) EnsureDir(JoinPath(sub, b.folder));

              // 3x3 table (shape optional)
              Make3x3Table_TH1(ds, histBase, sub,
                               string("table3x3_") + outStem + (shapeTable ? "_shape.png" : ".png"),
                               xTitle, shapeTable ? "A.U." : "Counts",
                               logy, shapeTable, common);

              // per-bin plots
              for (int i = 0; i < kNPtBins; ++i)
              {
                const PtBin& b = PtBins()[i];
                TH1* h = GetObj<TH1>(ds, histBase + b.suffix, true, true, true);
                if (!h) continue;

                TH1* hc = CloneTH1(h, TString::Format("%s_%s_%d", histBase.c_str(), outStem.c_str(), i).Data());
                if (!hc) continue;

                if (shapeTable) NormalizeToUnitArea(hc);

                vector<string> lines = common;
                lines.push_back(TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data());
                lines.push_back(TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                lines.push_back(TString::Format("Fiducial: |#eta_{jet}| < %.1f", etaFidAbs).Data());

                const string fp = JoinPath(sub, b.folder + "/" + outStem + "_" + b.folder + (shapeTable ? "_shape.png" : ".png"));
                DrawAndSaveTH1_Common(ds, hc, fp, xTitle, shapeTable ? "A.U." : "Counts", lines, logy, drawFidLines, etaFidAbs);

                delete hc;
              }
            };

            // --- jet1 family ---
            plotVar("h_jet1Pt_" + rKey,      dirJet1, "jet1Pt",   "p_{T}^{jet1} [GeV]", true,  false, false);
            plotVar("h_jet1Eta_sel_" + rKey, dirJet1, "jet1Eta",  "#eta_{jet1}",        false, false, true);
            plotVar("h_jet1Phi_sel_" + rKey, dirJet1, "jet1Phi",  "#phi_{jet1}",        false, false, false);
            plotVar("h_jet1Mass_sel_" + rKey,dirJet1, "jet1Mass", "m_{jet1} [GeV]",     false, false, false);

            // --- jet2 family ---
            plotVar("h_jet2Pt_" + rKey,      dirJet2, "jet2Pt",   "p_{T}^{jet2} [GeV]", true,  false, false);
            plotVar("h_jet2Eta_sel_" + rKey, dirJet2, "jet2Eta",  "#eta_{jet2}",        false, false, true);
            plotVar("h_jet2Phi_sel_" + rKey, dirJet2, "jet2Phi",  "#phi_{jet2}",        false, false, false);

            // Per-bin counts + mean pT summary (terminal + file)
            for (int i = 0; i < kNPtBins; ++i)
            {
              const PtBin& b = PtBins()[i];

              TH1* h1 = GetObj<TH1>(ds, string("h_jet1Pt_") + rKey + b.suffix, true, true, true);
              TH1* h2 = GetObj<TH1>(ds, string("h_jet2Pt_") + rKey + b.suffix, true, true, true);

              const double n1 = h1 ? h1->GetEntries() : 0.0;
              const double n2 = h2 ? h2->GetEntries() : 0.0;
              const double m1 = h1 ? h1->GetMean()    : 0.0;
              const double m2 = h2 ? h2->GetMean()    : 0.0;

              cout << std::left << std::setw(wPt) << b.label
                   << std::right
                   << std::setw(wN) << std::fixed << std::setprecision(0) << n1
                   << std::setw(wN) << n2
                   << std::setw(wM) << std::fixed << std::setprecision(3) << m1
                   << std::setw(wM) << m2
                   << "\n";

              summaryLines.push_back(TString::Format(
                "pT=%s  Njet1=%.0f  Njet2=%.0f  <pT1>=%.6f  <pT2>=%.6f",
                b.label.c_str(), n1, n2, m1, m2
              ).Data());
            }

            WriteTextFile(JoinPath(rOut, "summary_selectedJets.txt"), summaryLines);
          }
        }

        // =============================================================================
        // Section 5E: xJ + alpha QA (per pTgamma bin) + terminal summaries + overlays
        // Uses:
        //   h_xJ_rXX_pT_lo_hi
        //   h_alpha_rXX_pT_lo_hi
        // =============================================================================
        void RunXJAlphaQA(Dataset& ds)
        {
            cout << ANSI_BOLD_CYN << "\n==============================\n"
                 << "[SECTION 5E] x_{J} and #alpha QA (" << ds.label << ")\n"
                 << "==============================" << ANSI_RESET << "\n";

            string outDir = ds.isSim
              ? JoinPath(ds.outBase, "RecoilJetQA/xJAlpha")
              : JoinPath(ds.outBase, "baselineData/RecoilJetQA/xJAlpha");

            EnsureDir(outDir);

            vector<string> common;
            common.push_back("Selected recoil-jet kinematics (after matching)");
            common.push_back("Jets: p_{T}^{jet} #geq 10 GeV");

            // ---------------------------------------------------------------------------
            // Local helpers (refactor only: same plots, same filenames, same styles)
            // ---------------------------------------------------------------------------

            struct XJAlphaDirs
            {
              string rOut;
              string dirXJ;
              string dirAlpha;
              string dirOv;
            };

            auto PrintRKeySummaryHeader =
              [&](const string& rKey, double R)
            {
              cout << ANSI_BOLD_CYN
                   << "\n[xJ/alpha summary] " << ds.label << "  rKey=" << rKey
                   << " (R=" << std::fixed << std::setprecision(1) << R << ")\n"
                   << ANSI_RESET;
            };

            auto PrintTableHeader =
              [&](int wPt, int wN, int wM)
            {
              cout << std::left << std::setw(wPt) << "pTbin"
                   << std::right
                   << std::setw(wN) << "N"
                   << std::setw(wM) << "<xJ>"
                   << std::setw(wM) << "<alpha>"
                   << "\n";
              cout << string(wPt + wN + 2*wM, '-') << "\n";
            };

            auto InitSummaryLines =
              [&](const string& rKey, double R)->vector<string>
            {
              vector<string> lines;
              lines.push_back(string("xJ/alpha summary (") + ds.label + ")");
              lines.push_back(string("rKey: ") + rKey + TString::Format("  R=%.1f", R).Data());
              lines.push_back("");
              return lines;
            };

            auto MakeDirsForRKey =
              [&](const string& rKey)->XJAlphaDirs
            {
              XJAlphaDirs D;
              D.rOut     = JoinPath(outDir, rKey);
              D.dirXJ    = JoinPath(D.rOut, "xJ");
              D.dirAlpha = JoinPath(D.rOut, "alpha");
              D.dirOv    = JoinPath(D.rOut, "overlays");

              EnsureDir(D.rOut);
              EnsureDir(D.dirXJ);
              EnsureDir(D.dirAlpha);
              EnsureDir(D.dirOv);

              for (const auto& b : PtBins())
              {
                EnsureDir(JoinPath(D.dirXJ, b.folder));
                EnsureDir(JoinPath(D.dirAlpha, b.folder));
              }

              return D;
            };

            auto MakeShapeTablesForRKey =
              [&](const string& rKey, const XJAlphaDirs& D)
            {
              Make3x3Table_TH1(ds, "h_xJ_" + rKey, D.dirXJ,
                               "table3x3_xJ_shape.png",
                               "x_{J}", "A.U.",
                               false, true, common);

              Make3x3Table_TH1(ds, "h_alpha_" + rKey, D.dirAlpha,
                               "table3x3_alpha_shape.png",
                               "#alpha", "A.U.",
                               false, true, common);
            };

            struct MeansPack
            {
              vector<double> xCenters;
              vector<double> meanXJ;
              vector<double> meanA;
            };

            auto FillPerBinPlotsAndMeans =
              [&](const string& rKey, double R, const XJAlphaDirs& D,
                  int wPt, int wN, int wM,
                  vector<string>& summaryLines)->MeansPack
            {
              MeansPack M;
              M.xCenters.assign(kNPtBins, 0.0);
              M.meanXJ.assign(kNPtBins, 0.0);
              M.meanA.assign(kNPtBins, 0.0);

              for (int i = 0; i < kNPtBins; ++i)
              {
                const PtBin& b = PtBins()[i];
                M.xCenters[i] = 0.5*(b.lo + b.hi);

                TH1* hx = GetObj<TH1>(ds, "h_xJ_" + rKey + b.suffix, true, true, true);
                TH1* ha = GetObj<TH1>(ds, "h_alpha_" + rKey + b.suffix, true, true, true);

                const double N  = hx ? hx->GetEntries() : 0.0;
                const double mx = hx ? hx->GetMean()    : 0.0;
                const double ma = ha ? ha->GetMean()    : 0.0;

                M.meanXJ[i] = mx;
                M.meanA[i]  = ma;

                cout << std::left << std::setw(wPt) << b.label
                     << std::right
                     << std::setw(wN) << std::fixed << std::setprecision(0) << N
                     << std::setw(wM) << std::fixed << std::setprecision(4) << mx
                     << std::setw(wM) << std::fixed << std::setprecision(4) << ma
                     << "\n";

                summaryLines.push_back(TString::Format(
                  "pT=%s  N=%.0f  <xJ>=%.6f  <alpha>=%.6f",
                  b.label.c_str(), N, mx, ma
                ).Data());

                // xJ plot (shape)
                if (hx)
                {
                  TH1* hc = CloneTH1(hx, TString::Format("xJ_%s_%d", rKey.c_str(), i).Data());
                  if (hc)
                  {
                    NormalizeToUnitArea(hc);
                    vector<string> extra = common;
                    extra.push_back(TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data());
                    extra.push_back(TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                    DrawAndSaveTH1_Common(ds, hc,
                      JoinPath(D.dirXJ, b.folder + "/xJ_shape_" + b.folder + ".png"),
                      "x_{J}", "A.U.", extra, false);
                    delete hc;
                  }
                }

                // alpha plot (shape)
                if (ha)
                {
                  TH1* hc = CloneTH1(ha, TString::Format("alpha_%s_%d", rKey.c_str(), i).Data());
                  if (hc)
                  {
                    NormalizeToUnitArea(hc);
                    vector<string> extra = common;
                    extra.push_back(TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data());
                    extra.push_back(TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                    DrawAndSaveTH1_Common(ds, hc,
                      JoinPath(D.dirAlpha, b.folder + "/alpha_shape_" + b.folder + ".png"),
                      "#alpha", "A.U.", extra, false);
                    delete hc;
                  }
                }
              }

              return M;
            };

            auto DrawMeanOverlaysForRKey =
              [&](const string& rKey, const XJAlphaDirs& D, const MeansPack& M)
            {
              // mean xJ vs pTgamma
              {
                TCanvas c1(TString::Format("c_meanxJ_%s", rKey.c_str()).Data(), "c_meanxJ", 900, 700);
                ApplyCanvasMargins1D(c1);
                TGraph g(kNPtBins, &M.xCenters[0], &M.meanXJ[0]);
                g.SetLineWidth(2);
                g.SetMarkerStyle(20);
                g.Draw("ALP");
                g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g.GetYaxis()->SetTitle("<x_{J}>");
                g.SetMinimum(0.0);
                g.SetMaximum(1.2);

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {string("Mean x_{J} vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);
                SaveCanvas(c1, JoinPath(D.dirOv, "mean_xJ_vs_pTgamma.png"));
              }

              // mean alpha vs pTgamma
              {
                TCanvas c2(TString::Format("c_meana_%s", rKey.c_str()).Data(), "c_meana", 900, 700);
                ApplyCanvasMargins1D(c2);
                TGraph g(kNPtBins, &M.xCenters[0], &M.meanA[0]);
                g.SetLineWidth(2);
                g.SetMarkerStyle(20);
                g.Draw("ALP");
                g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g.GetYaxis()->SetTitle("<#alpha>");
                g.SetMinimum(0.0);
                g.SetMaximum(1.0);

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {string("Mean #alpha vs p_{T}^{#gamma}"), rKey}, 0.030, 0.040);
                SaveCanvas(c2, JoinPath(D.dirOv, "mean_alpha_vs_pTgamma.png"));
              }
            };

            auto ComputeMeansAcrossRKeys =
              [&](map<string, vector<double> >& meanXJ,
                  map<string, vector<double> >& meanA)
            {
              meanXJ.clear();
              meanA.clear();

              for (const auto& rKey : kRKeys)
              {
                meanXJ[rKey] = vector<double>(kNPtBins, 0.0);
                meanA[rKey]  = vector<double>(kNPtBins, 0.0);

                for (int i = 0; i < kNPtBins; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  TH1* hx = GetObj<TH1>(ds, "h_xJ_" + rKey + b.suffix, false, false, false);
                  TH1* ha = GetObj<TH1>(ds, "h_alpha_" + rKey + b.suffix, false, false, false);
                  meanXJ[rKey][i] = hx ? hx->GetMean() : 0.0;
                  meanA[rKey][i]  = ha ? ha->GetMean() : 0.0;
                }
              }
            };

            auto DrawRKeyOverlayGraphs =
              [&](const string& overDir,
                  const vector<double>& x,
                  const map<string, vector<double> >& meanXJ,
                  const map<string, vector<double> >& meanA)
            {
              // mean xJ overlay
              {
                TCanvas c("c_ov_meanxJ","c_ov_meanxJ",900,700);
                ApplyCanvasMargins1D(c);

                TGraph g02(kNPtBins, &x[0], &meanXJ.at("r02")[0]);
                TGraph g04(kNPtBins, &x[0], &meanXJ.at("r04")[0]);
                g02.SetLineWidth(2); g04.SetLineWidth(2);
                g02.SetMarkerStyle(20); g04.SetMarkerStyle(24);
                g02.SetLineColor(1); g04.SetLineColor(2);
                g02.SetMarkerColor(1); g04.SetMarkerColor(2);

                g02.Draw("ALP");
                g02.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g02.GetYaxis()->SetTitle("<x_{J}>");
                g02.SetMinimum(0.0);
                g02.SetMaximum(1.2);
                g04.Draw("LP same");

                TLegend leg(0.62,0.78,0.92,0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.033);
                leg.AddEntry(&g02, "r02 (R=0.2)", "lp");
                leg.AddEntry(&g04, "r04 (R=0.4)", "lp");
                leg.Draw();

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {"Overlay: <x_{J}> vs p_{T}^{#gamma}"}, 0.030, 0.040);

                SaveCanvas(c, JoinPath(overDir, "overlay_mean_xJ_r02_vs_r04.png"));
              }

              // mean alpha overlay
              {
                TCanvas c("c_ov_meana","c_ov_meana",900,700);
                ApplyCanvasMargins1D(c);

                TGraph g02(kNPtBins, &x[0], &meanA.at("r02")[0]);
                TGraph g04(kNPtBins, &x[0], &meanA.at("r04")[0]);
                g02.SetLineWidth(2); g04.SetLineWidth(2);
                g02.SetMarkerStyle(20); g04.SetMarkerStyle(24);
                g02.SetLineColor(1); g04.SetLineColor(2);
                g02.SetMarkerColor(1); g04.SetMarkerColor(2);

                g02.Draw("ALP");
                g02.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g02.GetYaxis()->SetTitle("<#alpha>");
                g02.SetMinimum(0.0);
                g02.SetMaximum(1.0);
                g04.Draw("LP same");

                TLegend leg(0.62,0.78,0.92,0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.033);
                leg.AddEntry(&g02, "r02 (R=0.2)", "lp");
                leg.AddEntry(&g04, "r04 (R=0.4)", "lp");
                leg.Draw();

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {"Overlay: <#alpha> vs p_{T}^{#gamma}"}, 0.030, 0.040);

                SaveCanvas(c, JoinPath(overDir, "overlay_mean_alpha_r02_vs_r04.png"));
              }
            };

            // ---------------------------------------------------------------------------
            // 1) Per-rKey: tables, per-bin plots, per-rKey overlays, text summaries
            // ---------------------------------------------------------------------------

            for (const auto& rKey : kRKeys)
            {
              const double R = RFromKey(rKey);

              PrintRKeySummaryHeader(rKey, R);

              const int wPt = 10;
              const int wN  = 12;
              const int wM  = 12;

              PrintTableHeader(wPt, wN, wM);

              vector<string> lines = InitSummaryLines(rKey, R);

              const XJAlphaDirs D = MakeDirsForRKey(rKey);

              MakeShapeTablesForRKey(rKey, D);

              MeansPack M = FillPerBinPlotsAndMeans(rKey, R, D, wPt, wN, wM, lines);

              WriteTextFile(JoinPath(D.rOut, "summary_xJ_alpha.txt"), lines);

              DrawMeanOverlaysForRKey(rKey, D, M);
            }

            // ---------------------------------------------------------------------------
            // 2) Overlay across rKeys (r02 vs r04) for mean xJ and mean alpha
            // ---------------------------------------------------------------------------

            {
              map<string, vector<double> > meanXJ;
              map<string, vector<double> > meanA;
              ComputeMeansAcrossRKeys(meanXJ, meanA);

              vector<double> x(kNPtBins, 0.0);
              for (int i = 0; i < kNPtBins; ++i) x[i] = 0.5*(PtBins()[i].lo + PtBins()[i].hi);

              const string overDir = JoinPath(outDir, "overlays");
              EnsureDir(overDir);

              DrawRKeyOverlayGraphs(overDir, x, meanXJ, meanA);
            }
        }
  
  
  
        // =============================================================================
        // JES3 In-situ residual calibration (DATA vs SIM) using reco-only JES3 TH3s.
        //
        // This function is designed to be called once per dataset from inside RunJES3QA.
        // It internally caches pointers to the SIM dataset and all DATA datasets seen,
        // and only runs the calibration once SIM+DATA are both available.
        //
        // Enabled ONLY when:
        //   inline bool isSimAndDataPP = true;
        //
        // Output + behavior matches the previously inlined lambda block.
        // =============================================================================
        void JES3_InSituResidualCalibration_MaybeRun(Dataset& ds)
        {
            if (!isSimAndDataPP) return;

            // ---- Cross-call cache (RunJES3QA is invoked once per dataset) ----
            struct InSituCache
            {
              Dataset* sim = nullptr;                   // chosen SIM dataset
              std::vector<Dataset*> data;              // all DATA datasets seen
              std::set<std::string> doneDataLabels;    // which DATA labels already processed
            };
            static InSituCache C;

            auto ContainsAny = [&](const std::string& s, const std::vector<std::string>& needles)->bool
            {
              for (const auto& n : needles)
              {
                if (s.find(n) != std::string::npos) return true;
              }
              return false;
            };

            auto IsPreferredMergedSIM = [&](const Dataset& s)->bool
            {
              const SimSample sel = CurrentSimSample();
              if (!IsMergedSimSample(sel)) return true;

              // Prefer the dataset whose path matches the *selected* merged sample path
              const std::string want = SimInputPathForSample(sel);
              if (!want.empty() && s.inFilePath == want) return true;

              const std::string hay =
                s.label + " " + s.topDirName + " " + s.inFilePath;

              const std::string need = SimSampleLabel(sel);
              if (!need.empty() && hay.find(need) != std::string::npos) return true;

              return ContainsAny(hay, {"Merged","merged","MERGED","5and10","5and20","10and20","5and10and20"});
            };

            auto RegisterDataset = [&](Dataset& cur)
            {
              if (cur.isSim)
              {
                if (!C.sim)
                {
                  C.sim = &cur;
                }
                else
                {
                  const bool oldPref = IsPreferredMergedSIM(*C.sim);
                  const bool newPref = IsPreferredMergedSIM(cur);
                  if (newPref && !oldPref)
                  {
                    C.sim = &cur;
                  }
                }
              }
              else
              {
                bool exists = false;
                for (auto* p : C.data)
                {
                  if (p && p->label == cur.label) { exists = true; break; }
                }
                if (!exists) C.data.push_back(&cur);
              }
            };

            RegisterDataset(ds);

            // Need both SIM and at least one DATA dataset before we can do anything.
            if (!C.sim) return;
            if (C.data.empty()) return;

            // ---- Peak fit helper (iterative Gaussian) ----
            struct PeakFitResult
            {
              bool   ok    = false;
              double mu    = 0.0;
              double sigma = 0.0;

              double x0    = 0.0;   // initial max-bin center guess

              double r1Lo  = 0.0;
              double r1Hi  = 0.0;
              double r2Lo  = 0.0;
              double r2Hi  = 0.0;

              double chi2  = 0.0;
              int    ndf   = 0;
            };

            auto ClampRangeToAxis = [&](const TAxis* ax, double& lo, double& hi)
            {
              if (!ax) return;
              const double xmin = ax->GetXmin();
              const double xmax = ax->GetXmax();
              if (lo < xmin) lo = xmin;
              if (hi > xmax) hi = xmax;
              if (hi <= lo)
              {
                lo = xmin;
                hi = xmax;
              }
            };

            auto FitPeakIterativeGaussian = [&](TH1* h)->PeakFitResult
            {
              PeakFitResult R;
              if (!h) return R;
              if (h->GetEntries() <= 0.0) return R;

              const TAxis* ax = h->GetXaxis();
              if (!ax) return R;

              const int bMax = h->GetMaximumBin();
              R.x0 = ax->GetBinCenter(bMax);

              // Fit #1: [x0-0.15, x0+0.15]
              R.r1Lo = R.x0 - 0.15;
              R.r1Hi = R.x0 + 0.15;
              ClampRangeToAxis(ax, R.r1Lo, R.r1Hi);

              TF1 f1("f1_gaus","gaus", R.r1Lo, R.r1Hi);
              f1.SetLineWidth(2);
              h->Fit(&f1, "RQ0S");

              const double mu1 = f1.GetParameter(1);
              const double s1  = std::fabs(f1.GetParameter(2));
              const double sForRange = (s1 > 0.0 ? s1 : 0.10);

              // Fit #2: [mu1-1.5*sigma1, mu1+1.5*sigma1]
              R.r2Lo = mu1 - 1.5*sForRange;
              R.r2Hi = mu1 + 1.5*sForRange;
              ClampRangeToAxis(ax, R.r2Lo, R.r2Hi);

              TF1 f2("f2_gaus","gaus", R.r2Lo, R.r2Hi);
              f2.SetLineWidth(2);
              h->Fit(&f2, "RQ0S");

              R.mu    = f2.GetParameter(1);
              R.sigma = std::fabs(f2.GetParameter(2));
              R.chi2  = f2.GetChisquare();
              R.ndf   = f2.GetNDF();

              R.ok = std::isfinite(R.mu) && std::isfinite(R.sigma) && (R.sigma > 0.0);
              return R;
            };

            // ---- Drawing helper: 1D hist with Gaussian overlay and annotation ----
            auto DrawXJWithFit =
              [&](Dataset& outDs,
                  TH1* h,
                  const PeakFitResult& fit,
                  const std::string& outPng,
                  const std::string& xTitle,
                  const std::vector<std::string>& lines,
                  int colorHist,
                  int colorFit)
            {
              if (!h) return;

              TCanvas c(TString::Format("c_insitu_%s", outPng.c_str()).Data(), "c_insitu", 900, 700);
              ApplyCanvasMargins1D(c);

              h->SetTitle("");
              h->SetLineWidth(2);
              h->SetMarkerStyle(20);
              h->SetMarkerSize(1.00);
              h->SetLineColor(colorHist);
              h->SetMarkerColor(colorHist);
              h->SetFillStyle(0);

              h->GetXaxis()->SetTitle(xTitle.c_str());
              h->GetYaxis()->SetTitle("Counts");

              const double ymax = h->GetMaximum();
              if (ymax > 0.0) h->SetMaximum(ymax * 1.25);

              h->Draw("E1");

              if (fit.ok)
              {
                TF1 f("f_gaus_draw","gaus", fit.r2Lo, fit.r2Hi);
                f.SetLineColor(colorFit);
                f.SetLineWidth(2);
                f.SetParameters(h->GetMaximum(), fit.mu, fit.sigma);
                f.Draw("same");

                TLine l1(fit.r2Lo, 0.0, fit.r2Lo, h->GetMaximum());
                TLine l2(fit.r2Hi, 0.0, fit.r2Hi, h->GetMaximum());
                l1.SetLineStyle(2); l2.SetLineStyle(2);
                l1.SetLineColor(colorFit); l2.SetLineColor(colorFit);
                l1.Draw("same"); l2.Draw("same");
              }

              DrawLatexLines(0.14, 0.92, DefaultHeaderLines(outDs), 0.034, 0.045);
              DrawLatexLines(0.14, 0.80, lines, 0.030, 0.040);

              SaveCanvas(c, outPng);
            };

            // ---- Drawing helper: overlay DATA vs MC (shape) ----
            auto DrawOverlayDataVsMCShape =
              [&](Dataset& outDs,
                  TH1* hData,
                  TH1* hMC,
                  const std::string& outPng,
                  const std::vector<std::string>& lines)
            {
              if (!hData || !hMC) return;

              TH1* d = CloneTH1(hData, "hData_shape_tmp");
              TH1* m = CloneTH1(hMC,   "hMC_shape_tmp");
              if (!d || !m)
              {
                if (d) delete d;
                if (m) delete m;
                return;
              }

              NormalizeToUnitArea(d);
              NormalizeToUnitArea(m);

              TCanvas c(TString::Format("c_ov_insitu_%s", outPng.c_str()).Data(), "c_ov_insitu", 900, 700);
              ApplyCanvasMargins1D(c);

              d->SetTitle("");
              d->SetLineWidth(2);
              d->SetMarkerStyle(20);
              d->SetMarkerSize(1.00);
              d->SetLineColor(1);
              d->SetMarkerColor(1);

              m->SetLineWidth(2);
              m->SetMarkerStyle(24);
              m->SetMarkerSize(1.00);
              m->SetLineColor(2);
              m->SetMarkerColor(2);

              d->GetXaxis()->SetTitle("x_{J#gamma}");
              d->GetYaxis()->SetTitle("A.U.");

              const double ymax = std::max(d->GetMaximum(), m->GetMaximum());
              if (ymax > 0.0) d->SetMaximum(ymax * 1.25);

              d->Draw("E1");
              m->Draw("E1 same");

              TLegend leg(0.62, 0.78, 0.92, 0.90);
              leg.SetTextFont(42);
              leg.SetTextSize(0.033);
              leg.AddEntry(d, "DATA (shape)", "ep");
              leg.AddEntry(m, "MC (shape)",   "ep");
              leg.Draw();

              DrawLatexLines(0.14, 0.92, DefaultHeaderLines(outDs), 0.034, 0.045);
              DrawLatexLines(0.14, 0.80, lines, 0.030, 0.040);

              SaveCanvas(c, outPng);

              delete d;
              delete m;
            };
        }
  
  
        void JES3_R02R04Overlays_MaybeRun(Dataset& ds, const std::string& outDir)
        {
              if (!ds.isSim) return;

              const std::string ovBase = JoinPath(outDir, "r02_r04_r06");
              EnsureDir(ovBase);

              auto AlphaTag = [&](double aMax)->std::string
              {
                  std::ostringstream s;
                  s << std::fixed << std::setprecision(2) << aMax;
                  std::string t = s.str();
                  std::replace(t.begin(), t.end(), '.', 'p');
                  return std::string("alphaLT") + t;
              };

              auto DrawOverlayPair_TH3xJ =
                [&](const TH3* h02, const TH3* h04, const TH3* h06,
                    const std::string& outDirHere,
                    const std::string& xTitle,
                    const std::vector<std::string>& headerLines,
                    bool useAlphaCut,
                    double alphaMax)
              {
                  if (!h02 || !h04) return;

                  EnsureDir(outDirHere);

                  const int n02 = h02->GetXaxis()->GetNbins();
                  const int n04 = h04->GetXaxis()->GetNbins();
                  int nPt = std::min(n02, n04);
                  if (h06)
                  {
                      const int n06 = h06->GetXaxis()->GetNbins();
                      nPt = std::min(nPt, n06);
                  }

                  // Per-bin overlay PNGs
                  for (int ib = 1; ib <= nPt; ++ib)
                  {
                      TH1* a = nullptr;
                      TH1* b = nullptr;
                      TH1* c6 = nullptr;

                      if (useAlphaCut)
                      {
                          a = ProjectY_AtXbin_AndAlphaMax_TH3(
                                h02, ib, alphaMax,
                                TString::Format("xJ_ov_r02_alphaLT%.2f_b%d", alphaMax, ib).Data());

                          b = ProjectY_AtXbin_AndAlphaMax_TH3(
                                h04, ib, alphaMax,
                                TString::Format("xJ_ov_r04_alphaLT%.2f_b%d", alphaMax, ib).Data());

                          if (h06)
                          {
                              c6 = ProjectY_AtXbin_AndAlphaMax_TH3(
                                     h06, ib, alphaMax,
                                     TString::Format("xJ_ov_r06_alphaLT%.2f_b%d", alphaMax, ib).Data());
                          }
                      }
                      else
                      {
                          a = ProjectY_AtXbin_TH3(
                                h02, ib,
                                TString::Format("xJ_ov_r02_int_b%d", ib).Data());

                          b = ProjectY_AtXbin_TH3(
                                h04, ib,
                                TString::Format("xJ_ov_r04_int_b%d", ib).Data());

                          if (h06)
                          {
                              c6 = ProjectY_AtXbin_TH3(
                                     h06, ib,
                                     TString::Format("xJ_ov_r06_int_b%d", ib).Data());
                          }
                      }

                      if (!a && !b && !c6)
                      {
                          if (a) delete a;
                          if (b) delete b;
                          if (c6) delete c6;
                          continue;
                      }

                      const std::string ptLab = AxisBinLabel(h02->GetXaxis(), ib, "GeV", 0);

                      TCanvas can(TString::Format("c_ov_%s_%d", outDirHere.c_str(), ib).Data(), "c_ov", 900, 700);
                      ApplyCanvasMargins1D(can);

                      if (a)
                      {
                          a->SetDirectory(nullptr);
                          EnsureSumw2(a);

                          a->SetTitle("");
                          a->SetLineWidth(2);
                          a->SetLineColor(2);
                          a->SetMarkerStyle(20);
                          a->SetMarkerSize(1.05);
                          a->SetMarkerColor(2);
                          a->SetFillStyle(0);

                          a->GetXaxis()->SetTitle(xTitle.c_str());
                          a->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");
                      }

                      if (b)
                      {
                          b->SetDirectory(nullptr);
                          EnsureSumw2(b);

                          b->SetTitle("");
                          b->SetLineWidth(2);
                          b->SetLineColor(4);
                          b->SetMarkerStyle(20);
                          b->SetMarkerSize(1.05);
                          b->SetMarkerColor(4);
                          b->SetFillStyle(0);

                          b->GetXaxis()->SetTitle(xTitle.c_str());
                          b->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");
                      }

                      if (c6)
                      {
                          c6->SetDirectory(nullptr);
                          EnsureSumw2(c6);

                          c6->SetTitle("");
                          c6->SetLineWidth(2);
                          c6->SetLineColor(8);
                          c6->SetMarkerStyle(20);
                          c6->SetMarkerSize(1.05);
                          c6->SetMarkerColor(8);
                          c6->SetFillStyle(0);

                          c6->GetXaxis()->SetTitle(xTitle.c_str());
                          c6->GetYaxis()->SetTitle((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts");
                      }

                      double ymax = 0.0;
                      if (a)  ymax = std::max(ymax, a->GetMaximum());
                      if (b)  ymax = std::max(ymax, b->GetMaximum());
                      if (c6) ymax = std::max(ymax, c6->GetMaximum());

                      TH1* first = a ? a : (b ? b : c6);
                      if (first) first->SetMaximum(ymax * 1.25);

                      if (a) a->Draw("E1");
                      else if (b) b->Draw("E1");
                      else if (c6) c6->Draw("E1");

                      if (b)  b->Draw("E1 same");
                      if (c6) c6->Draw("E1 same");

                      TLegend leg(0.70, 0.78, 0.92, 0.90);
                      leg.SetTextFont(42);
                      leg.SetTextSize(0.035);
                      if (a)  leg.AddEntry(a,  "R = 0.2", "ep");
                      if (b)  leg.AddEntry(b,  "R = 0.4", "ep");
                      if (c6) leg.AddEntry(c6, "R = 0.6", "ep");
                      leg.Draw();

                      DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);

                      std::vector<std::string> lines = headerLines;
                      lines.push_back(TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());
                      if (useAlphaCut) lines.push_back(TString::Format("#alpha < %.2f", alphaMax).Data());
                      DrawLatexLines(0.14, 0.875, lines, 0.030, 0.038);

                      const auto& cfg = DefaultSim10and20Config();
                      std::vector<std::string> cutLines;
                      cutLines.push_back(TString::Format("Back-to-back: #Delta#phi_{#gamma,jet} > %s", cfg.bbLabel.c_str()).Data());
                      cutLines.push_back(TString::Format("p_{T}^{jet} > %.0f GeV", cfg.jetMinPt).Data());
                      DrawLatexLines(0.70, 0.74, cutLines, 0.030, 0.038);

                      SaveCanvas(can, JoinPath(outDirHere, TString::Format("overlay_pTbin%d.png", ib).Data()));

                      if (a) delete a;
                      if (b) delete b;
                      if (c6) delete c6;
                  }

                  // Single 3x2 summary table using ONLY the last 6 JES3 pT bins (skip first bin)
                  const int nCols   = 3;
                  const int nRows   = 2;
                  const int perPage = nCols * nRows; // 6

                  const int startBinForTable = 2; // skip 13-15; start at 15-17
                  const int nTableBins = (nPt >= startBinForTable) ? std::min(perPage, nPt - startBinForTable + 1) : 0;

                  if (nTableBins > 0)
                  {
                      TCanvas canTbl(
                        TString::Format("c_tbl_%s_last6", outDirHere.c_str()).Data(),
                        "c_tbl_overlay_last6", 1500, 900
                      );
                      canTbl.Divide(nCols, nRows, 0.001, 0.001);

                      std::vector<TH1*> keep;
                      keep.reserve(3 * perPage);

                      for (int k = 0; k < nTableBins; ++k)
                      {
                          const int ib = startBinForTable + k;
                          canTbl.cd(k + 1);

                          gPad->SetLeftMargin(0.14);
                          gPad->SetRightMargin(0.05);
                          gPad->SetBottomMargin(0.14);
                          gPad->SetTopMargin(0.10);
                          gPad->SetGrid(0,0);

                          TH1* a = nullptr;
                          TH1* b = nullptr;
                          TH1* c6 = nullptr;

                          if (useAlphaCut)
                          {
                              a = ProjectY_AtXbin_AndAlphaMax_TH3(
                                    h02, ib, alphaMax,
                                    TString::Format("tbl_r02_alphaLT%.2f_b%d", alphaMax, ib).Data());

                              b = ProjectY_AtXbin_AndAlphaMax_TH3(
                                    h04, ib, alphaMax,
                                    TString::Format("tbl_r04_alphaLT%.2f_b%d", alphaMax, ib).Data());

                              if (h06)
                              {
                                  c6 = ProjectY_AtXbin_AndAlphaMax_TH3(
                                         h06, ib, alphaMax,
                                         TString::Format("tbl_r06_alphaLT%.2f_b%d", alphaMax, ib).Data());
                              }
                          }
                          else
                          {
                              a = ProjectY_AtXbin_TH3(
                                    h02, ib,
                                    TString::Format("tbl_r02_int_b%d", ib).Data());

                              b = ProjectY_AtXbin_TH3(
                                    h04, ib,
                                    TString::Format("tbl_r04_int_b%d", ib).Data());

                              if (h06)
                              {
                                  c6 = ProjectY_AtXbin_TH3(
                                         h06, ib,
                                         TString::Format("tbl_r06_int_b%d", ib).Data());
                              }
                          }

                          if (!a && !b && !c6)
                          {
                              if (a) delete a;
                              if (b) delete b;
                              if (c6) delete c6;
                              continue;
                          }

                          const bool doShapeNorm = (!useAlphaCut && (outDirHere.find("RECO") != std::string::npos));

                          if (a)
                          {
                              a->SetDirectory(nullptr);
                              EnsureSumw2(a);

                              a->SetTitle("");
                              a->SetLineWidth(2);
                              a->SetLineColor(2);
                              a->SetMarkerStyle(20);
                              a->SetMarkerSize(1.05);
                              a->SetMarkerColor(2);
                              a->SetFillStyle(0);

                              a->GetXaxis()->SetTitle(xTitle.c_str());
                              a->GetYaxis()->SetTitle(doShapeNorm ? "A.U." : ((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts"));
                          }

                          if (b)
                          {
                              b->SetDirectory(nullptr);
                              EnsureSumw2(b);

                              b->SetTitle("");
                              b->SetLineWidth(2);
                              b->SetLineColor(4);
                              b->SetMarkerStyle(20);
                              b->SetMarkerSize(1.05);
                              b->SetMarkerColor(4);
                              b->SetFillStyle(0);

                              b->GetXaxis()->SetTitle(xTitle.c_str());
                              b->GetYaxis()->SetTitle(doShapeNorm ? "A.U." : ((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts"));
                          }

                          if (c6)
                          {
                              c6->SetDirectory(nullptr);
                              EnsureSumw2(c6);

                              c6->SetTitle("");
                              c6->SetLineWidth(2);
                              c6->SetLineColor(8);
                              c6->SetMarkerStyle(20);
                              c6->SetMarkerSize(1.05);
                              c6->SetMarkerColor(8);
                              c6->SetFillStyle(0);

                              c6->GetXaxis()->SetTitle(xTitle.c_str());
                              c6->GetYaxis()->SetTitle(doShapeNorm ? "A.U." : ((ds.isSim && IsWeightedSIMSelected()) ? "Counts / pb^{-1}" : "Counts"));
                          }

                          if (doShapeNorm)
                          {
                              if (a)  NormalizeToUnitArea(a);
                              if (b)  NormalizeToUnitArea(b);
                              if (c6) NormalizeToUnitArea(c6);
                          }

                          const double xMaxPlot = 2.0;
                          if (a)  a->GetXaxis()->SetRangeUser(0.0, xMaxPlot);
                          if (b)  b->GetXaxis()->SetRangeUser(0.0, xMaxPlot);
                          if (c6) c6->GetXaxis()->SetRangeUser(0.0, xMaxPlot);

                          double ymax = 0.0;
                          if (a)  ymax = std::max(ymax, a->GetMaximum());
                          if (b)  ymax = std::max(ymax, b->GetMaximum());
                          if (c6) ymax = std::max(ymax, c6->GetMaximum());

                          TH1* first = a ? a : (b ? b : c6);
                          if (first) first->SetMaximum(ymax * 1.08);

                          if (a) a->Draw("E1");
                          else if (b) b->Draw("E1");
                          else if (c6) c6->Draw("E1");

                          if (b)  b->Draw("E1 same");
                          if (c6) c6->Draw("E1 same");

                          const std::string ptLab = AxisBinLabel(h02->GetXaxis(), ib, "GeV", 0);

                          TLatex ttitle;
                          ttitle.SetNDC(true);
                          ttitle.SetTextFont(42);
                          ttitle.SetTextAlign(22);
                          ttitle.SetTextSize(0.060);

                          const bool isTruthPlot = (std::string(xTitle).find("truth") != std::string::npos);
                          const char* levelTag   = isTruthPlot ? "Truth-level" : "Reco-level";
                          const char* pTTag      = isTruthPlot ? "p_{T}^{#gamma,truth}" : "p_{T}^{#gamma}";

                          ttitle.DrawLatex(
                            0.50, 0.95,
                            TString::Format("%s %s, %s = %s", levelTag, xTitle.c_str(), pTTag, ptLab.c_str()).Data()
                          );

                          const auto& cfgDef = DefaultSim10and20Config();

                          TLatex tcuts;
                          tcuts.SetNDC(true);
                          tcuts.SetTextFont(42);
                          tcuts.SetTextAlign(33);
                          tcuts.SetTextSize(0.040);
                          tcuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", cfgDef.bbLabel.c_str()).Data());
                          tcuts.DrawLatex(0.92, 0.54, TString::Format("Reco (p_{T}^{jet} > %.0f GeV)", cfgDef.jetMinPt).Data());

                          TLegend leg(0.72, 0.68, 0.92, 0.90);
                          leg.SetTextFont(42);
                          leg.SetTextSize(0.065);
                          leg.SetFillStyle(0);
                          leg.SetBorderSize(0);
                          if (a)  leg.AddEntry(a,  "R = 0.2", "ep");
                          if (b)  leg.AddEntry(b,  "R = 0.4", "ep");
                          if (c6) leg.AddEntry(c6, "R = 0.6", "ep");
                          leg.DrawClone();

                          if (useAlphaCut)
                          {
                              TLatex ta;
                              ta.SetNDC(true);
                              ta.SetTextFont(42);
                              ta.SetTextAlign(13);
                              ta.SetTextSize(0.050);
                              ta.DrawLatex(0.80, 0.70, TString::Format("#alpha < %.2f", alphaMax).Data());
                          }

                          if (a)  keep.push_back(a);
                          if (b)  keep.push_back(b);
                          if (c6) keep.push_back(c6);
                      }

                      const std::string outName = useAlphaCut
                        ? TString::Format("table3x2_overlay_alphaLT%.2f.png", alphaMax).Data()
                        : "table3x2_overlay_integratedAlpha.png";

                      SaveCanvas(canTbl, JoinPath(outDirHere, outName));

                      for (auto* h : keep) delete h;
                  }
              };

              const std::string dirInt  = JoinPath(ovBase, "xJ_integratedAlpha");
              const std::string dirCuts = JoinPath(ovBase, "xJ_alphaCuts");
              EnsureDir(dirInt);
              EnsureDir(dirCuts);

              // -------------------- TRUTH (reco-conditioned, jet-matched): integrated alpha --------------------
              TH3* hTr02 = GetObj<TH3>(ds, "h_JES3Truth_pT_xJ_alpha_r02", true, true, true);
              TH3* hTr04 = GetObj<TH3>(ds, "h_JES3Truth_pT_xJ_alpha_r04", true, true, true);

              if (hTr02 && hTr04)
              {
                  const std::string outHere = JoinPath(dirInt, "TRUTH_recoConditioned");
                  DrawOverlayPair_TH3xJ(
                    hTr02, hTr04, nullptr,
                    outHere,
                    "x_{J#gamma}^{truth}",
                    {"TRUTH (reco-conditioned, jet-matched): x_{J#gamma}^{truth}", "Overlay: r02 red, r04 blue"},
                    false, 0.0
                  );
              }
              else
              {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] TRUTH reco-conditioned overlay skipped: missing h_JES3Truth_pT_xJ_alpha_r02 or r04 in dataset " << ds.label
                       << ANSI_RESET << "\n";
              }

              // -------------------- TRUTH (pure): integrated alpha --------------------
              TH3* hTrPure02 = GetObj<TH3>(ds, "h_JES3TruthPure_pT_xJ_alpha_r02", true, true, true);
              TH3* hTrPure04 = GetObj<TH3>(ds, "h_JES3TruthPure_pT_xJ_alpha_r04", true, true, true);

              if (hTrPure02 && hTrPure04)
              {
                  const std::string outHere = JoinPath(dirInt, "TRUTH_pure");
                  DrawOverlayPair_TH3xJ(
                    hTrPure02, hTrPure04, nullptr,
                    outHere,
                    "x_{J#gamma}^{truth}",
                    {"TRUTH (pure): x_{J#gamma}^{truth}", "Overlay: r02 red, r04 blue"},
                    false, 0.0
                  );
              }
              else
              {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] TRUTH pure overlay skipped: missing h_JES3TruthPure_pT_xJ_alpha_r02 or r04 in dataset " << ds.label
                       << ANSI_RESET << "\n";
              }

              // -------------------- RECO (baseline): integrated alpha + alpha cuts --------------------
              TH3* hRe02 = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_r02", true, true, true);
              TH3* hRe04 = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_r04", true, true, true);
              TH3* hRe06 = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_r06", true, true, true);

              if (hRe02 && hRe04 && hRe06)
              {
                  const std::string outHere = JoinPath(dirInt, "RECO");
                  DrawOverlayPair_TH3xJ(
                    hRe02, hRe04, hRe06,
                    outHere,
                    "x_{J#gamma}",
                    {"RECO: x_{J#gamma}", "Overlay: r02 red, r04 blue, r06 dark green"},
                    false, 0.0
                  );

                  const std::vector<double> alphaMaxCuts = {0.20, 0.30, 0.40, 0.50};
                  const std::string cutsBase = JoinPath(dirCuts, "RECO");
                  EnsureDir(cutsBase);

                  for (double aMax : alphaMaxCuts)
                  {
                      const std::string aDir = JoinPath(cutsBase, AlphaTag(aMax));
                      DrawOverlayPair_TH3xJ(
                        hRe02, hRe04, hRe06,
                        aDir,
                        "x_{J#gamma}",
                        {"RECO: x_{J#gamma}", "Overlay: r02 red, r04 blue, r06 dark green"},
                        true, aMax
                      );
                  }
              }
              else
              {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] RECO overlays skipped: missing h_JES3_pT_xJ_alpha_r02 or r04 or r06 in dataset " << ds.label
                       << ANSI_RESET << "\n";
              }

              // -------------------- RECO (truth-PHOTON-tagged): integrated alpha --------------------
              TH3* hRePhoTag02 = GetObj<TH3>(ds, "h_JES3RecoTruthPhoTagged_pT_xJ_alpha_r02", true, true, true);
              TH3* hRePhoTag04 = GetObj<TH3>(ds, "h_JES3RecoTruthPhoTagged_pT_xJ_alpha_r04", true, true, true);

              if (hRePhoTag02 && hRePhoTag04)
              {
                  const std::string outHere = JoinPath(dirInt, "RECO_truthPhoTagged");
                  DrawOverlayPair_TH3xJ(
                    hRePhoTag02, hRePhoTag04, nullptr,
                    outHere,
                    "x_{J#gamma}",
                    {"RECO (truth-PHOTON-tagged): x_{J#gamma}", "Overlay: r02 red, r04 blue"},
                    false, 0.0
                  );
              }
              else
              {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] RECO truth-PHOTON-tagged overlay skipped: missing h_JES3RecoTruthPhoTagged_pT_xJ_alpha_r02 or r04 in dataset " << ds.label
                       << ANSI_RESET << "\n";
              }

              // -------------------- RECO (truth-tagged PHOTON+JET): integrated alpha --------------------
              TH3* hReTag02 = GetObj<TH3>(ds, "h_JES3RecoTruthTagged_pT_xJ_alpha_r02", true, true, true);
              TH3* hReTag04 = GetObj<TH3>(ds, "h_JES3RecoTruthTagged_pT_xJ_alpha_r04", true, true, true);

              if (hReTag02 && hReTag04)
              {
                  const std::string outHere = JoinPath(dirInt, "RECO_truthTaggedPhoJet");
                  DrawOverlayPair_TH3xJ(
                    hReTag02, hReTag04, nullptr,
                    outHere,
                    "x_{J#gamma}",
                    {"RECO (truth-tagged PHOTON+JET): x_{J#gamma}", "Overlay: r02 red, r04 blue"},
                    false, 0.0
                  );
              }
              else
              {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] RECO truth-tagged (PHOTON+JET) overlay skipped: missing h_JES3RecoTruthTagged_pT_xJ_alpha_r02 or r04 in dataset " << ds.label
                       << ANSI_RESET << "\n";
              }

              // -------------------- TRUTH (reco-conditioned, NO jet match): integrated alpha --------------------
              TH3* hTrNoJM02 = GetObj<TH3>(ds, "h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_r02", true, true, true);
              TH3* hTrNoJM04 = GetObj<TH3>(ds, "h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_r04", true, true, true);

              if (hTrNoJM02 && hTrNoJM04)
              {
                  const std::string outHere = JoinPath(dirInt, "TRUTH_recoConditioned_noJetMatch");
                  DrawOverlayPair_TH3xJ(
                    hTrNoJM02, hTrNoJM04, nullptr,
                    outHere,
                    "x_{J#gamma}^{truth}",
                    {"TRUTH (reco-conditioned, NO jet match): x_{J#gamma}^{truth}", "Overlay: r02 red, r04 blue"},
                    false, 0.0
                  );
              }
              else
              {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] TRUTH reco-conditioned (NO jet match) overlay skipped: missing h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_r02 or r04 in dataset " << ds.label
                       << ANSI_RESET << "\n";
              }

              // =============================================================================
              // recoSampleOverlays (RECO JES3 xJ overlays across samples)
              //
              // Runs ONLY when the selected SIM sample is:
              //   allPhoton5and10and20sim  ->  SimSample::kPhotonJet5And10And20Merged
              //
              // Output folder structure:
              //   <outDir>/recoSampleOverlays/r02/
              //   <outDir>/recoSampleOverlays/r04/
              //
              // For EACH radius folder, produces THREE PNGs (3x3 table across pT bins):
              //   1) 7-curve overlay: 5,10,20, (5+10), (5+20), (10+20), (5+10+20)
              //   2) 4-curve overlay: 5,10,20, (5+10+20)
              //   3) 3-curve overlay: (5+10), (10+20), (5+10+20)
              // =============================================================================
              if (ds.isSim && CurrentSimSample() == SimSample::kPhotonJet5And10And20Merged)
              {
                  cout << ANSI_BOLD_CYN
                       << "\n[JES3 recoSampleOverlays] Building RECO xJ overlay tables across photonJet samples...\n"
                       << ANSI_RESET;

                  auto EnsureMerged = [&](const std::string& outFile,
                                          const std::vector<std::string>& ins,
                                          const std::vector<double>& sigmas,
                                          const std::vector<std::string>& labs)->bool
                  {
                      if (!gSystem->AccessPathName(outFile.c_str())) return true;
                      return BuildMergedSIMFile_PhotonSlices(ins, sigmas, outFile, kDirSIM, labs);
                  };

                  EnsureMerged(kMergedSIMOut_5and10,
                               {kInSIM5, DefaultSim10and20Config().photon10},
                               {kSigmaPhoton5_pb, kSigmaPhoton10_pb},
                               {"photonJet5", "photonJet10"});

                  EnsureMerged(kMergedSIMOut_5and20,
                               {kInSIM5, DefaultSim10and20Config().photon20},
                               {kSigmaPhoton5_pb, kSigmaPhoton20_pb},
                               {"photonJet5", "photonJet20"});

                  EnsureMerged(MergedSIMOut_10and20_Default(),
                               {DefaultSim10and20Config().photon10, DefaultSim10and20Config().photon20},
                               {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                               {"photonJet10", "photonJet20"});

                  const std::string base = JoinPath(outDir, "recoSampleOverlays");
                  EnsureDir(base);

                  struct Sample
                  {
                      std::string legend;
                      std::string filePath;
                      int color = 1;
                      bool isMerged = false;
                      TFile* f = nullptr;
                      TDirectory* top = nullptr;
                  };

                  std::vector<Sample> S =
                  {
                      {"photon jet 5 GeV",            kInSIM5,                             kBlack,      false, nullptr, nullptr},
                      {"photon jet 10 GeV",           DefaultSim10and20Config().photon10,  kRed+1,      false, nullptr, nullptr},
                      {"photon jet 20 GeV",           DefaultSim10and20Config().photon20,  kBlue+1,     false, nullptr, nullptr},
                      {"photon jet 5 + 10 GeV",       kMergedSIMOut_5and10,                kGreen+3,    true,  nullptr, nullptr},
                      {"photon jet 5 + 20 GeV",       kMergedSIMOut_5and20,                kMagenta+1,  true,  nullptr, nullptr},
                      {"photon jet 10 + 20 GeV",      MergedSIMOut_10and20_Default(),      kOrange+7,   true,  nullptr, nullptr},
                      {"photon jet 5 + 10 + 20 GeV",  kMergedSIMOut_5and10and20,           kViolet+1,   true,  nullptr, nullptr}
                  };

                  auto CloseAll = [&]()
                  {
                      for (auto& s : S)
                      {
                          if (s.f)
                          {
                              s.f->Close();
                              s.f = nullptr;
                              s.top = nullptr;
                          }
                      }
                  };

                  auto OpenSample = [&](Sample& s)->bool
                  {
                      s.f = TFile::Open(s.filePath.c_str(), "READ");
                      if (!s.f || s.f->IsZombie())
                      {
                          cout << ANSI_BOLD_YEL
                               << "[WARN] recoSampleOverlays: cannot open: " << s.filePath
                               << ANSI_RESET << "\n";
                          if (s.f) { s.f->Close(); s.f = nullptr; }
                          s.top = nullptr;
                          return false;
                      }

                      s.top = s.f->GetDirectory(kDirSIM.c_str());
                      if (!s.top)
                      {
                          cout << ANSI_BOLD_YEL
                               << "[WARN] recoSampleOverlays: missing topDir '" << kDirSIM << "' in: " << s.filePath
                               << ANSI_RESET << "\n";
                          s.f->Close();
                          s.f = nullptr;
                          s.top = nullptr;
                          return false;
                      }
                      return true;
                  };

                  for (auto& s : S) OpenSample(s);

                  auto GetTH3FromSample = [&](const Sample& s, const std::string& hname)->TH3*
                  {
                      if (!s.top) return nullptr;
                      return dynamic_cast<TH3*>(s.top->Get(hname.c_str()));
                  };

                  auto StyleOverlayHist =
                    [&](TH1* h,
                        const Sample& samp,
                        bool forceClosedCircles,
                        bool forceOpenCircles,
                        int  colorOverride)
                  {
                      if (!h) return;
                      h->SetDirectory(nullptr);
                      EnsureSumw2(h);

                      const int nb = h->GetNbinsX();
                      const double area = h->Integral(1, nb, "width");
                      if (area > 0.0) h->Scale(1.0 / area);

                      const int col = (colorOverride >= 0 ? colorOverride : samp.color);

                      int mStyle = samp.isMerged ? 24 : 20;
                      if (forceClosedCircles) mStyle = 20;
                      if (forceOpenCircles)   mStyle = 24;

                      h->SetTitle("");
                      h->SetLineWidth(3);
                      h->SetLineColor(col);

                      h->SetMarkerStyle(mStyle);
                      h->SetMarkerSize(1.05);
                      h->SetMarkerColor(col);

                      h->SetFillStyle(0);

                      h->GetXaxis()->SetTitle("x_{J#gamma}");
                      h->GetYaxis()->SetTitle("Normalized (area = 1)");

                      h->GetYaxis()->SetTitleSize(0.075);
                      h->GetYaxis()->SetLabelSize(0.060);
                      h->GetYaxis()->SetTitleOffset(0.95);

                      h->GetXaxis()->SetTitleSize(0.070);
                      h->GetXaxis()->SetLabelSize(0.055);
                  };

                  auto Make3x3OverlayTable =
                    [&](const std::string& rKey,
                        const std::string& outDirR,
                        const std::string& outStem,
                        const std::vector<int>& useIdx,
                        bool forceClosedCircles = false,
                        bool forceOpenCircles   = false)
                  {
                      EnsureDir(outDirR);

                      const std::string th3Name = "h_JES3_pT_xJ_alpha_" + rKey;

                      std::vector<TH3*> H3;
                      H3.reserve(useIdx.size());

                      const TAxis* axPt = nullptr;
                      for (int idx : useIdx)
                      {
                          TH3* h3 = GetTH3FromSample(S[idx], th3Name);
                          H3.push_back(h3);
                          if (!axPt && h3) axPt = h3->GetXaxis();
                      }

                      if (!axPt)
                      {
                          cout << ANSI_BOLD_YEL
                               << "[WARN] recoSampleOverlays: missing " << th3Name << " for all requested samples (rKey=" << rKey << ")"
                               << ANSI_RESET << "\n";
                          return;
                      }

                      const int nPt = axPt->GetNbins();
                      const int perPage = 9;
                      int page = 0;

                      for (int start = 1; start <= nPt; start += perPage)
                      {
                          ++page;

                          TCanvas can(
                            TString::Format("c_recoSampleOv_%s_%s_p%d", rKey.c_str(), outStem.c_str(), page).Data(),
                            "c_recoSampleOv", 1500, 1200
                          );
                          can.Divide(3, 3, 0.001, 0.001);

                          std::vector<TH1*> keep;
                          keep.reserve(perPage * useIdx.size());

                          for (int k = 0; k < perPage; ++k)
                          {
                              const int ib = start + k;
                              can.cd(k + 1);

                              gPad->SetLeftMargin(0.14);
                              gPad->SetRightMargin(0.05);
                              gPad->SetBottomMargin(0.14);
                              gPad->SetTopMargin(0.10);
                              gPad->SetLogy(false);

                              if (ib > nPt)
                              {
                                  TLatex t;
                                  t.SetNDC(true);
                                  t.SetTextFont(42);
                                  t.SetTextSize(0.06);
                                  t.DrawLatex(0.20, 0.55, "EMPTY");
                                  continue;
                              }

                              const std::string ptLab = AxisBinLabel(axPt, ib, "GeV", 0);

                              std::vector<TH1*> proj;
                              proj.reserve(useIdx.size());

                              double ymax = 0.0;

                              for (size_t j = 0; j < useIdx.size(); ++j)
                              {
                                  TH3* h3 = H3[j];
                                  if (!h3)
                                  {
                                      proj.push_back(nullptr);
                                      continue;
                                  }

                                  const int sIdx = useIdx[j];

                                  TH1* hx = ProjectY_AtXbin_TH3(
                                    h3, ib,
                                    TString::Format("xJ_recoSampleOv_%s_%s_b%d_s%d", rKey.c_str(), outStem.c_str(), ib, sIdx).Data()
                                  );

                                  if (!hx)
                                  {
                                      proj.push_back(nullptr);
                                      continue;
                                  }

                                  int colorOverride = -1;
                                  if (outStem == "recoJES3_xJ_overlays_10and20_vs_5and10and20")
                                  {
                                      if (sIdx == 5) colorOverride = kBlue+1;
                                      if (sIdx == 6) colorOverride = kRed+1;
                                  }

                                  StyleOverlayHist(hx, S[sIdx], forceClosedCircles, forceOpenCircles, colorOverride);
                                  ymax = std::max(ymax, hx->GetMaximum());

                                  proj.push_back(hx);
                                  keep.push_back(hx);
                              }

                              TH1* first = nullptr;
                              for (TH1* h : proj) { if (h) { first = h; break; } }

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
                              for (TH1* h : proj)
                              {
                                  if (!h) continue;
                                  if (!drawn) { h->Draw("E1"); drawn = true; }
                                  else        { h->Draw("E1 same"); }
                              }

                              TLatex t;
                              t.SetNDC(true);
                              t.SetTextFont(42);
                              t.SetTextAlign(13);
                              t.SetTextSize(0.062);
                              t.DrawLatex(0.14, 0.965, TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());

                              const bool isTwoCurveSpecial =
                                (outStem == "recoJES3_xJ_overlays_10and20_vs_5and10and20");

                              TLegend leg(
                                isTwoCurveSpecial ? 0.38 : 0.48,
                                isTwoCurveSpecial ? 0.73 : 0.52,
                                isTwoCurveSpecial ? 0.76 : 0.93,
                                isTwoCurveSpecial ? 0.88 : 0.92
                              );

                              leg.SetTextFont(42);
                              leg.SetTextSize(isTwoCurveSpecial ? 0.054 : (useIdx.size() >= 6 ? 0.042 : 0.052));
                              leg.SetFillStyle(0);
                              leg.SetBorderSize(0);
                              leg.SetMargin(0.25);

                              for (size_t j = 0; j < proj.size(); ++j)
                              {
                                  if (!proj[j]) continue;
                                  leg.AddEntry(proj[j], S[useIdx[j]].legend.c_str(), "ep");
                              }
                              leg.DrawClone();
                          }

                          const std::string outName =
                            (nPt <= perPage)
                              ? TString::Format("table3x3_%s.png", outStem.c_str()).Data()
                              : TString::Format("table3x3_%s_page%d.png", outStem.c_str(), page).Data();

                          SaveCanvas(can, JoinPath(outDirR, outName));

                          for (auto* h : keep) delete h;
                      }
                  };

                  for (const auto& rKey : kRKeys)
                  {
                      const std::string outR = JoinPath(base, rKey);
                      EnsureDir(outR);

                      Make3x3OverlayTable(
                        rKey, outR,
                        "recoJES3_xJ_overlays_allSamples",
                        {0, 1, 2, 3, 4, 5, 6}
                      );

                      Make3x3OverlayTable(
                        rKey, outR,
                        "recoJES3_xJ_overlays_singlesOnly",
                        {0, 1, 2},
                        true,
                        false
                      );

                      Make3x3OverlayTable(
                        rKey, outR,
                        "recoJES3_xJ_overlays_singlesPlusAllMerged",
                        {0, 1, 2, 6}
                      );

                      Make3x3OverlayTable(
                        rKey, outR,
                        "recoJES3_xJ_overlays_mergedPairsPlusAllMerged",
                        {3, 5, 6}
                      );

                      Make3x3OverlayTable(
                        rKey, outR,
                        "recoJES3_xJ_overlays_10and20_vs_5and10and20",
                        {5, 6},
                        false,
                        true
                      );
                  }

                  CloseAll();

                  cout << ANSI_BOLD_GRN
                       << "[OK] JES3 recoSampleOverlays written under: " << base << "\n"
                       << ANSI_RESET;
              }
        }

        // =============================================================================
        // JES3 RECO-only xJ overlay: compare Δφ cuts (π/2 vs 7π/8), r04 only
        //
        // Runs ONLY when:
        //   - ds.isSim == true
        //   - CurrentSimSample() == SimSample::kPhotonJet10And20Merged   (bothPhoton10and20sim=true)
        //
        // Baseline analysis continues using kInSIM10/kInSIM20 (pihalves).
        // This helper additionally reads the *_7piOver8.root slice files and builds an
        // in-memory weighted (10+20) merged TH3 for the single histogram:
        //   h_JES3_pT_xJ_alpha_r04
        //
        // Output:
        //   <outDir>/r04/xJ_fromJES3/RECO/table3x2_overlay_integratedAlpha_dPhiCuts.png
        //
        // Style:
        //   π/2    -> red, closed circles
        //   7π/8   -> blue, open circles
        // =============================================================================
        void JES3_RecoDeltaPhiCutOverlay_MaybeRun(Dataset& ds, const std::string& outDir)
        {
          if (!ds.isSim) return;
          if (CurrentSimSample() != SimSample::kPhotonJet10And20Merged) return;

          // We only want r04 for this request
          const std::string rKey = "r04";
          const std::string outRecoDir = JoinPath(JoinPath(outDir, rKey), "xJ_fromJES3/RECO");
          EnsureDir(outRecoDir);

          const std::string outPng =
            JoinPath(outRecoDir, "table3x2_overlay_integratedAlpha_dPhiCuts.png");

          // Baseline (π/2) comes from the currently-open merged dataset file
          TH3* hBase = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_r04", true, true, true);
          if (!hBase)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Δφ overlay skipped: missing h_JES3_pT_xJ_alpha_r04 in baseline merged dataset.\n"
                 << ANSI_RESET;
            return;
          }

            // Alternate (7π/8): MUST use the merged photonJet10+20 ROOT file for this cfgKey
            // IMPORTANT: read cfgKey from the header default so the jetMinPt printed + used is always consistent.
            const std::string altKey = DefaultSimSampleKey();
            const Sim10and20Config& altCfg = Sim10and20ConfigForKey(altKey);
            const std::string altMerged = MergedSIMOut_10and20_ForKey(altKey);

            const double jetMinPtGeV = altCfg.jetMinPt;
            const std::string bbLabel = altCfg.bbLabel;

            auto EnsureAltMerged = [&]()->bool
            {
              if (!gSystem->AccessPathName(altMerged.c_str())) return true; // exists

              cout << ANSI_BOLD_CYN
                   << "\n[MERGE] Building merged SIM10+20 file for Δφ overlay:\n"
                   << "  cfgKey   = " << altKey << "\n"
                   << "  in10     = " << altCfg.photon10 << "\n"
                   << "  in20     = " << altCfg.photon20 << "\n"
                   << "  out      = " << altMerged << "\n"
                   << ANSI_RESET;

              return BuildMergedSIMFile_PhotonSlices(
                {altCfg.photon10, altCfg.photon20},
                {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                altMerged,
                kDirSIM,
                {"photonJet10", "photonJet20"}
              );
            };

            if (!EnsureAltMerged())
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Δφ overlay skipped: could not build merged file: " << altMerged << "\n"
                   << ANSI_RESET;
              return;
            }

            TFile* fAlt = TFile::Open(altMerged.c_str(), "READ");
            if (!fAlt || fAlt->IsZombie())
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Δφ overlay skipped: cannot open merged file:\n"
                   << "  " << altMerged << "\n"
                   << ANSI_RESET;
              if (fAlt) fAlt->Close();
              return;
            }

            TDirectory* dAlt = fAlt->GetDirectory(kDirSIM.c_str());
            if (!dAlt)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Δφ overlay skipped: missing topDir '" << kDirSIM << "' in merged file:\n"
                   << "  " << altMerged << "\n"
                   << ANSI_RESET;
              fAlt->Close();
              return;
            }

            TH3* hAltIn = dynamic_cast<TH3*>(dAlt->Get("h_JES3_pT_xJ_alpha_r04"));
            if (!hAltIn)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Δφ overlay skipped: missing h_JES3_pT_xJ_alpha_r04 in merged file:\n"
                   << "  " << altMerged << "\n"
                   << ANSI_RESET;
              fAlt->Close();
              return;
            }

            TH3* hAlt = CloneTH3(hAltIn, "h_JES3_pT_xJ_alpha_r04_ALT_7piOver8");
            if (!hAlt)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Δφ overlay skipped: failed to clone alternate TH3 from merged file.\n"
                   << ANSI_RESET;
              fAlt->Close();
              return;
            }

            hAlt->SetDirectory(nullptr);
            if (hAlt->GetSumw2N() == 0) hAlt->Sumw2();

            fAlt->Close();

            cout << ANSI_DIM
                 << "  [Δφ overlay inputs]\n"
                 << "    baseline (opened dataset) = " << ds.inFilePath << "\n"
                 << "    alternate merged          = " << altMerged << "\n"
                 << ANSI_RESET;

            // -------------------------------------------------------------------------
            // Make a 3x2 table (6 pads) across pTgamma bins: ProjectY (xJ) integrating alpha
            // AND also write per-pT individual overlays (same styling) into perPtBin_dPhiCuts/
            // -------------------------------------------------------------------------
            const int nCols = 3;
            const int nRows = 2;

            const double R = RFromKey(rKey);

            const int nBase = hBase->GetXaxis()->GetNbins();
            const int nAlt  = hAlt->GetXaxis()->GetNbins();
            const int nPt   = std::min(nBase, nAlt);

            const std::string outPerPtDir = JoinPath(outRecoDir, "perPtBin_dPhiCuts");
            EnsureDir(outPerPtDir);

            TCanvas c("c_tbl_dphiCuts_r04", "c_tbl_dphiCuts_r04", 1500, 900);
            c.Divide(nCols, nRows, 0.001, 0.001);

            std::vector<TH1*> keep;
            keep.reserve(2 * 6);

            const int ibStart = 2; // skip first pT bin (e.g. 13-15) in the table
            const int nAvail  = (nPt >= ibStart) ? (nPt - ibStart + 1) : 0;
            const int nPads   = std::min(6, nAvail);

            for (int k = 0; k < nPads; ++k)
            {
              const int ib = ibStart + k;
              c.cd(k+1);

              gPad->SetLeftMargin(0.14);
              gPad->SetRightMargin(0.05);
              gPad->SetBottomMargin(0.14);
              gPad->SetTopMargin(0.10);
              gPad->SetLogy(false);

              TH1* hPi2  = ProjectY_AtXbin_TH3(hBase, ib, TString::Format("h_xJ_pi2_b%d", ib).Data());
              TH1* h7p8  = ProjectY_AtXbin_TH3(hAlt,  ib, TString::Format("h_xJ_7p8_b%d", ib).Data());

              if (!hPi2 && !h7p8)
              {
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(0.06);
                t.DrawLatex(0.15, 0.55, "MISSING");
                continue;
              }

                // Red closed circles: |Δφ| > π/2  (baseline)
                if (hPi2)
                {
                  hPi2->SetDirectory(nullptr);
                  EnsureSumw2(hPi2);
                  hPi2->SetTitle("");
                  hPi2->SetLineWidth(2);
                  hPi2->SetLineColor(2);
                  hPi2->SetMarkerStyle(20);   // closed circle
                  hPi2->SetMarkerSize(0.95);
                  hPi2->SetMarkerColor(2);
                  hPi2->SetFillStyle(0);
                  hPi2->GetXaxis()->SetTitle("x_{J#gamma}");
                  hPi2->GetYaxis()->SetTitle("A.U.");
                  hPi2->GetXaxis()->SetRangeUser(0.0, 2.0);
                  NormalizeToUnitArea(hPi2);
                }

                // Blue open circles: |Δφ| > 7π/8 (alternate)
                if (h7p8)
                {
                  h7p8->SetDirectory(nullptr);
                  EnsureSumw2(h7p8);
                  h7p8->SetTitle("");
                  h7p8->SetLineWidth(2);
                  h7p8->SetLineColor(4);
                  h7p8->SetMarkerStyle(24);   // open circle
                  h7p8->SetMarkerSize(0.95);
                  h7p8->SetMarkerColor(4);
                  h7p8->SetFillStyle(0);
                  h7p8->GetXaxis()->SetTitle("x_{J#gamma}");
                  h7p8->GetYaxis()->SetTitle("A.U.");
                  h7p8->GetXaxis()->SetRangeUser(0.0, 2.0);
                  NormalizeToUnitArea(h7p8);

                }

              TH1* first  = hPi2 ? hPi2 : h7p8;
              TH1* second = (first == hPi2) ? h7p8 : hPi2;

              double ymax = 0.0;
              if (hPi2) ymax = std::max(ymax, hPi2->GetMaximum());
              if (h7p8) ymax = std::max(ymax, h7p8->GetMaximum());
              if (first) first->SetMaximum(ymax * 1.25);

              if (first)  first->Draw("E1");
              if (second) second->Draw("E1 same");

              // Pad title: pT bin label taken from axis bin edges + include radius
              const std::string ptLab = AxisBinLabel(hBase->GetXaxis(), ib, "GeV", 0);

              TLatex ttitle;
              ttitle.SetNDC(true);
              ttitle.SetTextFont(42);
              ttitle.SetTextAlign(22);
              ttitle.SetTextSize(0.060);
              ttitle.DrawLatex(
                0.50, 0.95,
                TString::Format("Reco-level x_{J#gamma}, p_{T}^{#gamma} = %s  (R=%.1f)", ptLab.c_str(), R).Data()
              );

                // Legend: smaller and moved to the top-right (avoid marker overlap)
                TLegend leg(0.62, 0.74, 0.93, 0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.045);
                leg.SetFillStyle(0);
                leg.SetBorderSize(0);
                if (hPi2) leg.AddEntry(hPi2,  "|#Delta#phi(#gamma,jet)| > #pi/2",   "ep");
                if (h7p8) leg.AddEntry(h7p8,  "|#Delta#phi(#gamma,jet)| > 7#pi/8", "ep");
                leg.DrawClone();

                // Print the jetMinPt (from header default cfgKey) + back-to-back label under the legend
                TLatex tCuts;
                tCuts.SetNDC(true);
                tCuts.SetTextFont(42);
                tCuts.SetTextAlign(13); // left-top
                tCuts.SetTextSize(0.040);
                tCuts.DrawLatex(0.62, 0.71, TString::Format("p_{T}^{jet,min} = %.0f GeV", jetMinPtGeV).Data());
                tCuts.DrawLatex(0.62, 0.66, TString::Format("Back-to-back: %s", bbLabel.c_str()).Data());


              if (hPi2) keep.push_back(hPi2);
              if (h7p8) keep.push_back(h7p8);

              // ---------------------------------------------------------------------
              // ALSO write an individual per-pT overlay PNG with the same legend/header
              // ---------------------------------------------------------------------
              {
                TH1* pPi2 = ProjectY_AtXbin_TH3(hBase, ib, TString::Format("h_xJ_pi2_ind_b%d", ib).Data());
                TH1* p7p8 = ProjectY_AtXbin_TH3(hAlt,  ib, TString::Format("h_xJ_7p8_ind_b%d", ib).Data());

                if (pPi2 || p7p8)
                {
                    if (pPi2)
                    {
                      pPi2->SetDirectory(nullptr);
                      EnsureSumw2(pPi2);
                      pPi2->SetTitle("");
                      pPi2->SetLineWidth(2);
                      pPi2->SetLineColor(2);
                      pPi2->SetMarkerStyle(20);
                      pPi2->SetMarkerSize(1.00);
                      pPi2->SetMarkerColor(2);
                      pPi2->SetFillStyle(0);
                      pPi2->GetXaxis()->SetTitle("x_{J#gamma}");
                      pPi2->GetYaxis()->SetTitle("A.U.");
                      pPi2->GetXaxis()->SetRangeUser(0.0, 2.0);
                      NormalizeToUnitArea(pPi2);

                    }
                    if (p7p8)
                    {
                      p7p8->SetDirectory(nullptr);
                      EnsureSumw2(p7p8);
                      p7p8->SetTitle("");
                      p7p8->SetLineWidth(2);
                      p7p8->SetLineColor(4);
                      p7p8->SetMarkerStyle(24);
                      p7p8->SetMarkerSize(1.00);
                      p7p8->SetMarkerColor(4);
                      p7p8->SetFillStyle(0);
                      p7p8->GetXaxis()->SetTitle("x_{J#gamma}");
                      p7p8->GetYaxis()->SetTitle("A.U.");
                      p7p8->GetXaxis()->SetRangeUser(0.0, 2.0);
                      NormalizeToUnitArea(p7p8);
                    }

                  TH1* f1 = pPi2 ? pPi2 : p7p8;
                  TH1* f2 = (f1 == pPi2) ? p7p8 : pPi2;

                  double yM = 0.0;
                  if (pPi2) yM = std::max(yM, pPi2->GetMaximum());
                  if (p7p8) yM = std::max(yM, p7p8->GetMaximum());
                  if (f1) f1->SetMaximum(yM * 1.25);

                  TCanvas cInd(TString::Format("c_ind_dphiCuts_%s_b%d", rKey.c_str(), ib).Data(),
                              "c_ind_dphiCuts", 900, 700);
                  ApplyCanvasMargins1D(cInd);
                  cInd.SetLogy(false);

                  if (f1) f1->Draw("E1");
                  if (f2) f2->Draw("E1 same");

                  // Clean top-left header (same style as your perPt overlays)
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextAlign(13); // left-top

                  t.SetTextSize(0.036);
                  t.DrawLatex(0.14, 0.905, TString::Format("RECO: |#Delta#phi| cut comparison (%s)", rKey.c_str()).Data());

                  const std::string ptNoUnit = AxisBinLabel(hBase->GetXaxis(), ib, "", 0);
                  t.SetTextSize(0.032);
                  t.DrawLatex(0.14, 0.845,
                    TString::Format("p_{T}^{#gamma}: %s GeV (R = %.1f)", ptNoUnit.c_str(), R).Data()
                  );

                    // Legend (smaller, top-right to avoid marker overlap)
                    TLegend lInd(0.62, 0.78, 0.93, 0.92);
                    lInd.SetTextFont(42);
                    lInd.SetTextSize(0.032);
                    lInd.SetFillStyle(0);
                    lInd.SetBorderSize(0);
                    if (pPi2) lInd.AddEntry(pPi2, "|#Delta#phi(#gamma,jet)| > #pi/2",   "ep");
                    if (p7p8) lInd.AddEntry(p7p8, "|#Delta#phi(#gamma,jet)| > 7#pi/8", "ep");
                    lInd.Draw();

                    // Print jetMinPt + back-to-back label under the legend (from header default cfgKey)
                    t.SetTextSize(0.030);
                    t.DrawLatex(0.62, 0.74, TString::Format("p_{T}^{jet,min} = %.0f GeV", jetMinPtGeV).Data());
                    t.DrawLatex(0.62, 0.70, TString::Format("Back-to-back: %s", bbLabel.c_str()).Data());

                  const std::string outInd =
                    JoinPath(outPerPtDir, TString::Format("overlay_pTbin%d.png", ib).Data());
                  SaveCanvas(cInd, outInd);

                  if (pPi2) delete pPi2;
                  if (p7p8) delete p7p8;
                }
                else
                {
                  if (pPi2) delete pPi2;
                  if (p7p8) delete p7p8;
                }
              }
            }

            SaveCanvas(c, outPng);

            for (auto* h : keep) delete h;
            delete hAlt;

            cout << ANSI_BOLD_GRN
                 << "[OK] Wrote Δφ-cut RECO xJ overlay table: " << outPng << "\n"
                 << "[OK] Wrote Δφ-cut RECO per-pT overlays: " << outPerPtDir << "\n"
                 << ANSI_RESET;
        }

        // =============================================================================
        // JES3 RECO-only xJ overlay: compare p_{T}^{jet,min} (3 vs 5 vs 10), r02/r04/r06
        //
        // Runs ONLY when:
        //   - ds.isSim == true
        //   - CurrentSimSample() == SimSample::kPhotonJet10And20Merged   (bothPhoton10and20sim=true)
        //
        // Baseline (pTjet>10) comes from the currently-open merged dataset file:
        //   h_JES3_pT_xJ_alpha_<rKey>
        //
        // Additional curves are built in-memory by merging the corresponding slice files
        // with weights w = sigma / Naccepted, using the SAME back-to-back selection as
        // the current DefaultSimSampleKey() (pi/2 vs 7pi/8).
        //
        // Output (per radius):
        //   <outDir>/<rKey>/xJ_fromJES3/RECO/table3x2_overlay_integratedAlpha_pTminCompare.png
        //   <outDir>/<rKey>/xJ_fromJES3/RECO/perPtBin_pTminCompare/overlay_pTbinX.png
        //
        // Style (closed circles):
        //   pTjet > 3  -> black
        //   pTjet > 5  -> red
        //   pTjet > 10 -> blue
        // =============================================================================
        void JES3_RecoBinningOverlay_MaybeRun(Dataset& ds, const std::string& outDir)
        {
            if (!ds.isSim) return;
            if (CurrentSimSample() != SimSample::kPhotonJet10And20Merged) return;

            const bool is7pi8 = (DefaultSimSampleKey().find("7piOver8") != std::string::npos);

            const std::string keyPt5  = is7pi8 ? kAltSimSampleKey_jetMinPt5_7piOver8
                                               : kAltSimSampleKey_jetMinPt5_pihalves;
            const std::string keyPt3  = is7pi8 ? kAltSimSampleKey_jetMinPt3_7piOver8
                                               : kAltSimSampleKey_jetMinPt3_pihalves;

            const std::string bbLabel = Sim10and20ConfigForKey(keyPt3).bbLabel;

            auto BuildWeightedMergedTH3 =
              [&](const std::string& cfgKey, const std::string& h3name, const std::string& tag)->TH3*
            {
              const Sim10and20Config& cfg = Sim10and20ConfigForKey(cfgKey);
              const std::string merged = MergedSIMOut_10and20_ForKey(cfgKey);

              auto EnsureMerged = [&]()->bool
              {
                if (!gSystem->AccessPathName(merged.c_str())) return true; // exists

                cout << ANSI_BOLD_CYN
                     << "\n[MERGE] Building merged SIM10+20 file for pTmin overlay:\n"
                     << "  cfgKey   = " << cfgKey << "\n"
                     << "  in10     = " << cfg.photon10 << "\n"
                     << "  in20     = " << cfg.photon20 << "\n"
                     << "  out      = " << merged << "\n"
                     << ANSI_RESET;

                return BuildMergedSIMFile_PhotonSlices(
                  {cfg.photon10, cfg.photon20},
                  {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                  merged,
                  kDirSIM,
                  {"photonJet10", "photonJet20"}
                );
              };

              if (!EnsureMerged())
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] pTmin overlay skipped (" << tag << "): could not build merged file for cfgKey=" << cfgKey << "\n"
                     << "  out=" << merged << "\n"
                     << ANSI_RESET;
                return nullptr;
              }

              TFile* f = TFile::Open(merged.c_str(), "READ");
              if (!f || f->IsZombie())
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] pTmin overlay skipped (" << tag << "): cannot open merged file for cfgKey=" << cfgKey << "\n"
                     << "  " << merged << "\n"
                     << ANSI_RESET;
                if (f) f->Close();
                return nullptr;
              }

              TDirectory* d = f->GetDirectory(kDirSIM.c_str());
              if (!d)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] pTmin overlay skipped (" << tag << "): missing topDir '" << kDirSIM << "' in merged file for cfgKey=" << cfgKey << "\n"
                     << "  " << merged << "\n"
                     << ANSI_RESET;
                f->Close();
                return nullptr;
              }

              TH3* hin = dynamic_cast<TH3*>(d->Get(h3name.c_str()));
              if (!hin)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] pTmin overlay skipped (" << tag << "): missing " << h3name << " in merged file for cfgKey=" << cfgKey << "\n"
                     << "  " << merged << "\n"
                     << ANSI_RESET;
                f->Close();
                return nullptr;
              }

              TH3* h = CloneTH3(hin, TString::Format("%s_%s", h3name.c_str(), tag.c_str()).Data());
              if (h)
              {
                h->SetDirectory(nullptr);
                if (h->GetSumw2N() == 0) h->Sumw2();
              }

              f->Close();

              cout << ANSI_DIM
                   << "  [pTmin overlay input] cfgKey=" << cfgKey << "  merged=" << merged << "\n"
                   << ANSI_RESET;

              return h;
            };

            for (const auto& rKey : kRKeys)
            {
              const double R = RFromKey(rKey);

              const std::string outRecoDir = JoinPath(JoinPath(outDir, rKey), "xJ_fromJES3/RECO");
              EnsureDir(outRecoDir);

              const std::string outPng =
                JoinPath(outRecoDir, "table3x2_overlay_integratedAlpha_pTminCompare.png");

              const std::string h3name = "h_JES3_pT_xJ_alpha_" + rKey;

              // Baseline (pTjet > 10): comes from the currently-open merged dataset file
              TH3* hPt10 = GetObj<TH3>(ds, h3name, true, true, true);
              if (!hPt10)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] pTmin overlay skipped: missing " << h3name << " in baseline merged dataset.\n"
                     << ANSI_RESET;
                continue;
              }

              TH3* hPt5 = BuildWeightedMergedTH3(keyPt5, h3name, "pt5");
              TH3* hPt3 = BuildWeightedMergedTH3(keyPt3, h3name, "pt3");

              if (!hPt5 || !hPt3)
              {
                if (hPt5) delete hPt5;
                if (hPt3) delete hPt3;
                continue;
              }

              const int nCols = 3;
              const int nRows = 2;

              const int n10 = hPt10->GetXaxis()->GetNbins();
              const int n5  = hPt5->GetXaxis()->GetNbins();
              const int n3  = hPt3->GetXaxis()->GetNbins();
              const int nPt = std::min(n10, std::min(n5, n3));

              const std::string outPerPtDir = JoinPath(outRecoDir, "perPtBin_pTminCompare");
              EnsureDir(outPerPtDir);

              TCanvas c(TString::Format("c_tbl_pTminCompare_%s", rKey.c_str()).Data(),
                        TString::Format("c_tbl_pTminCompare_%s", rKey.c_str()).Data(),
                        1500, 900);
              c.Divide(nCols, nRows, 0.001, 0.001);

              std::vector<TH1*> keep;
              keep.reserve(3 * 6);

                auto Style = [&](TH1* h, int col)->void
                {
                  if (!h) return;
                  h->SetDirectory(nullptr);
                  EnsureSumw2(h);
                  h->SetTitle("");
                  h->SetLineWidth(2);
                  h->SetLineColor(col);
                  h->SetMarkerStyle(20);
                  h->SetMarkerSize(0.95);
                  h->SetMarkerColor(col);
                  h->SetFillStyle(0);
                  h->GetXaxis()->SetTitle("x_{J#gamma}");
                  h->GetXaxis()->SetRangeUser(0.0, 2.0);
                  h->GetYaxis()->SetTitle("A.U.");
                  NormalizeToUnitArea(h);
                };

                const int ibStart = 1; // start from first pT bin in the vector (ROOT bins are 1-indexed)
                const int nAvail  = (nPt >= ibStart) ? (nPt - ibStart + 1) : 0;
                const int nPads   = std::min(6, nAvail);

              for (int k = 0; k < nPads; ++k)
              {
                const int ib = ibStart + k;
                c.cd(k+1);

                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetLogy(false);

                TH1* h3XJ  = ProjectY_AtXbin_TH3(hPt3,  ib, TString::Format("h_xJ_pt3_%s_b%d",  rKey.c_str(), ib).Data());
                TH1* h5XJ  = ProjectY_AtXbin_TH3(hPt5,  ib, TString::Format("h_xJ_pt5_%s_b%d",  rKey.c_str(), ib).Data());
                TH1* h10XJ = ProjectY_AtXbin_TH3(hPt10, ib, TString::Format("h_xJ_pt10_%s_b%d", rKey.c_str(), ib).Data());

                if (!h3XJ && !h5XJ && !h10XJ)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.06);
                  t.DrawLatex(0.15, 0.55, "MISSING");
                  continue;
                }

                Style(h3XJ,  1); // black
                Style(h5XJ,  2); // red
                Style(h10XJ, 4); // blue

                TH1* base = h10XJ ? h10XJ : (h5XJ ? h5XJ : h3XJ);

                double ymax = 0.0;
                if (h3XJ)  ymax = std::max(ymax, h3XJ->GetMaximum());
                if (h5XJ)  ymax = std::max(ymax, h5XJ->GetMaximum());
                if (h10XJ) ymax = std::max(ymax, h10XJ->GetMaximum());
                if (base) base->SetMaximum(ymax * 1.25);

                if (base) base->Draw("E1");
                if (h3XJ  && h3XJ  != base) h3XJ->Draw("E1 same");
                if (h5XJ  && h5XJ  != base) h5XJ->Draw("E1 same");
                if (h10XJ && h10XJ != base) h10XJ->Draw("E1 same");

                const std::string ptLab = AxisBinLabel(hPt10->GetXaxis(), ib, "GeV", 0);

                TLatex ttitle;
                ttitle.SetNDC(true);
                ttitle.SetTextFont(42);
                ttitle.SetTextAlign(22);
                ttitle.SetTextSize(0.060);
                ttitle.DrawLatex(
                  0.50, 0.95,
                  TString::Format("Reco-level x_{J#gamma}, p_{T}^{#gamma} = %s  (R=%.1f)", ptLab.c_str(), R).Data()
                );

                TLatex tbb;
                tbb.SetNDC(true);
                tbb.SetTextFont(42);
                tbb.SetTextAlign(13);
                tbb.SetTextSize(0.050);
                tbb.DrawLatex(0.16, 0.86,
                  TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data()
                );

                TLegend leg(0.56, 0.62, 0.94, 0.86);
                leg.SetTextFont(42);
                leg.SetTextSize(0.040);
                leg.SetFillStyle(0);
                leg.SetBorderSize(0);
                if (h3XJ)  leg.AddEntry(h3XJ,  "Reco (p_{T}^{jet} > 3 GeV)",  "ep");
                if (h5XJ)  leg.AddEntry(h5XJ,  "Reco (p_{T}^{jet} > 5 GeV)",  "ep");
                if (h10XJ) leg.AddEntry(h10XJ, "Reco (p_{T}^{jet} > 10 GeV)", "ep");
                leg.DrawClone();

                if (h3XJ)  keep.push_back(h3XJ);
                if (h5XJ)  keep.push_back(h5XJ);
                if (h10XJ) keep.push_back(h10XJ);

                // Per-pT PNG
                {
                  TH1* p3  = ProjectY_AtXbin_TH3(hPt3,  ib, TString::Format("h_xJ_pt3_ind_%s_b%d",  rKey.c_str(), ib).Data());
                  TH1* p5  = ProjectY_AtXbin_TH3(hPt5,  ib, TString::Format("h_xJ_pt5_ind_%s_b%d",  rKey.c_str(), ib).Data());
                  TH1* p10 = ProjectY_AtXbin_TH3(hPt10, ib, TString::Format("h_xJ_pt10_ind_%s_b%d", rKey.c_str(), ib).Data());

                  if (p3 || p5 || p10)
                  {
                    Style(p3,  1);
                    Style(p5,  2);
                    Style(p10, 4);

                    TH1* f1 = p10 ? p10 : (p5 ? p5 : p3);

                    double yM = 0.0;
                    if (p3)  yM = std::max(yM, p3->GetMaximum());
                    if (p5)  yM = std::max(yM, p5->GetMaximum());
                    if (p10) yM = std::max(yM, p10->GetMaximum());
                    if (f1) f1->SetMaximum(yM * 1.25);

                    TCanvas cInd(TString::Format("c_ind_pTminCompare_%s_b%d", rKey.c_str(), ib).Data(),
                                "c_ind_pTminCompare", 900, 700);
                    ApplyCanvasMargins1D(cInd);
                    cInd.SetLogy(false);

                    if (f1) f1->Draw("E1");
                    if (p3  && p3  != f1) p3->Draw("E1 same");
                    if (p5  && p5  != f1) p5->Draw("E1 same");
                    if (p10 && p10 != f1) p10->Draw("E1 same");

                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(13);

                    t.SetTextSize(0.036);
                    t.DrawLatex(0.14, 0.905, TString::Format("RECO: p_{T}^{jet} min compare (%s)", rKey.c_str()).Data());

                    const std::string ptNoUnit = AxisBinLabel(hPt10->GetXaxis(), ib, "", 0);
                    t.SetTextSize(0.032);
                    t.DrawLatex(0.14, 0.845,
                      TString::Format("p_{T}^{#gamma}: %s GeV (R = %.1f)", ptNoUnit.c_str(), R).Data()
                    );

                    t.SetTextSize(0.032);
                    t.DrawLatex(0.14, 0.805,
                      TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data()
                    );

                    TLegend lInd(0.55, 0.72, 0.88, 0.90);
                    lInd.SetTextFont(42);
                    lInd.SetTextSize(0.040);
                    lInd.SetFillStyle(0);
                    lInd.SetBorderSize(0);
                    if (p3)  lInd.AddEntry(p3,  "p_{T}^{jet} > 3",  "ep");
                    if (p5)  lInd.AddEntry(p5,  "p_{T}^{jet} > 5",  "ep");
                    if (p10) lInd.AddEntry(p10, "p_{T}^{jet} > 10", "ep");
                    lInd.Draw();

                    const std::string outInd =
                      JoinPath(outPerPtDir, TString::Format("overlay_pTbin%d.png", ib).Data());
                    SaveCanvas(cInd, outInd);

                    if (p3)  delete p3;
                    if (p5)  delete p5;
                    if (p10) delete p10;
                  }
                  else
                  {
                    if (p3)  delete p3;
                    if (p5)  delete p5;
                    if (p10) delete p10;
                  }
                }
              }

              SaveCanvas(c, outPng);

              for (auto* h : keep) delete h;

              delete hPt5;
              delete hPt3;

              cout << ANSI_BOLD_GRN
                   << "[OK] Wrote pTmin RECO xJ overlay table: " << outPng << "\n"
                   << "[OK] Wrote pTmin RECO per-pT overlays: " << outPerPtDir << "\n"
                   << ANSI_RESET;
            }
      }

      void JES3_SamVsJustinUnsmearOverlay_MaybeRun(Dataset& ds, const std::string& outDir)
      {
          if (!ds.isSim) return;
          if (!doSamVsJustinUnsmearOverlays) return;
          (void)outDir; // this block is intentionally hard-coded to write ONLY under InputFilesSim/.../plots

          // Justin reference: jetMinPt3_7piOver8 (matched cuts)
          const Sim10and20Config& cfgJ = Sim10and20ConfigForKey(kAltSimSampleKey_jetMinPt3_7piOver8);

          const std::string baseDir  = DirFromPathSimple(cfgJ.photon10);
          const std::string plotsDir = JoinPath(baseDir, "plots");

          // Output: ONE folder only
          const std::string dirOver = JoinPath(plotsDir, "Overlays_RECO_SamVsJustin");

          // Two subfolders: normalized and unnormalized
          const std::string dirOver_norm = JoinPath(dirOver, "Normalized");
          const std::string dirOver_raw  = JoinPath(dirOver, "Unnormalized");

          // Sam inputs (ONLY these two files)
          const std::string sam10 = JoinPath(baseDir, "histsPhoton10.root");
          const std::string sam20 = JoinPath(baseDir, "histsPhoton20.root");

          EnsureDir(plotsDir);
          EnsureDir(dirOver);
          EnsureDir(dirOver_norm);
          EnsureDir(dirOver_raw);

          cout << ANSI_BOLD_CYN
               << "\n[SAM VS JUSTIN] Unsmear overlay (RECO-only, integrated over #alpha)\n"
               << "  baseDir    = " << baseDir << "\n"
               << "  plotsDir   = " << plotsDir << "\n"
               << "  overlays   = " << dirOver << "\n"
               << "  sam10      = " << sam10 << "\n"
               << "  sam20      = " << sam20 << "\n"
               << "  justin10   = " << cfgJ.photon10 << "\n"
               << "  justin20   = " << cfgJ.photon20 << "\n"
               << "  jetMinPt   = " << cfgJ.jetMinPt << " GeV\n"
               << "  backToBack = " << cfgJ.bbLabel << "\n"
               << ANSI_RESET;

          auto ReadXsecAndNev = [&](TFile* f, double& xsec_pb, double& nev_acc)->bool
          {
            xsec_pb = 0.0;
            nev_acc = 0.0;
            if (!f) return false;

            TH1* hx = dynamic_cast<TH1*>(f->Get("h_xsec"));
            if (!hx) hx = dynamic_cast<TH1*>(f->Get("xsec"));
            if (hx && hx->GetNbinsX() >= 1) xsec_pb = hx->GetBinContent(1);

            TH1* he = dynamic_cast<TH1*>(f->Get("h_evt_accept"));
            if (!he) he = dynamic_cast<TH1*>(f->Get("h_evtAccepted"));
            if (!he) he = dynamic_cast<TH1*>(f->Get("h_nevt"));
            if (!he) he = dynamic_cast<TH1*>(f->Get("h_evt"));
            if (he && he->GetNbinsX() >= 1) nev_acc = he->GetBinContent(1);

            return (xsec_pb > 0.0 && nev_acc > 0.0);
          };

          // -------------------------------------------------------------------------
          // Open Sam inputs (ONLY the two files below), and compute merge weights.
          // -------------------------------------------------------------------------
          TFile* fSam10 = TFile::Open(sam10.c_str(), "READ");
          TFile* fSam20 = TFile::Open(sam20.c_str(), "READ");
          if (!fSam10 || fSam10->IsZombie() || !fSam20 || fSam20->IsZombie())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Sam-vs-Justin overlay skipped: cannot open Sam files.\n"
                 << "       sam10=" << sam10 << "\n"
                 << "       sam20=" << sam20 << "\n"
                 << ANSI_RESET;
            if (fSam10) fSam10->Close();
            if (fSam20) fSam20->Close();
            return;
          }

          double sx10=0.0, sn10=0.0, sx20=0.0, sn20=0.0;
          const bool sok10 = ReadXsecAndNev(fSam10, sx10, sn10);
          const bool sok20 = ReadXsecAndNev(fSam20, sx20, sn20);

          double sw10 = 0.5;
          double sw20 = 0.5;
          if (sok10 && sok20)
          {
            sw10 = sx10 / sn10;
            sw20 = sx20 / sn20;
          }
          else
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] [SAM] Missing (xsec,nev) in Sam inputs; falling back to equal weights (0.5,0.5).\n"
                 << "       Expected histograms: h_xsec and h_evt_accept (bin1).\n"
                 << ANSI_RESET;
          }

          cout << ANSI_DIM
               << "[SAM] weights: sw10=" << std::setprecision(12) << sw10
               << "  sw20=" << std::setprecision(12) << sw20
               << "  (ok10=" << (sok10 ? "true" : "false") << ", ok20=" << (sok20 ? "true" : "false") << ")\n"
               << ANSI_RESET;

          // -------------------------------------------------------------------------
          // Open Justin slice files (jetMinPt3_7piOver8), and compute merge weights.
          // -------------------------------------------------------------------------
          TFile* fJ10 = TFile::Open(cfgJ.photon10.c_str(), "READ");
          TFile* fJ20 = TFile::Open(cfgJ.photon20.c_str(), "READ");
          if (!fJ10 || fJ10->IsZombie() || !fJ20 || fJ20->IsZombie())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Sam-vs-Justin overlay skipped: cannot open Justin slice files.\n"
                 << ANSI_RESET;
            if (fJ10) fJ10->Close();
            if (fJ20) fJ20->Close();
            fSam10->Close();
            fSam20->Close();
            return;
          }

          TDirectory* dJ10 = fJ10->GetDirectory(kDirSIM.c_str());
          TDirectory* dJ20 = fJ20->GetDirectory(kDirSIM.c_str());
          if (!dJ10 || !dJ20)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Sam-vs-Justin overlay skipped: missing topDir '" << kDirSIM << "' in Justin slice file(s).\n"
                 << ANSI_RESET;
            fJ10->Close();
            fJ20->Close();
            fSam10->Close();
            fSam20->Close();
            return;
          }

          const double N10 = ReadEventCountFromFile(fJ10, kDirSIM);
          const double N20 = ReadEventCountFromFile(fJ20, kDirSIM);
          if (N10 <= 0.0 || N20 <= 0.0)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Sam-vs-Justin overlay skipped: Naccepted <= 0 in Justin slice file(s).\n"
                 << ANSI_RESET;
            fJ10->Close();
            fJ20->Close();
            fSam10->Close();
            fSam20->Close();
            return;
          }

          const double jw10 = kSigmaPhoton10_pb / N10;
          const double jw20 = kSigmaPhoton20_pb / N20;

          cout << ANSI_DIM
               << "[JUSTIN] weights: jw10=" << std::setprecision(12) << jw10
               << "  jw20=" << std::setprecision(12) << jw20
               << "  (N10=" << std::fixed << std::setprecision(0) << N10
               << ", N20=" << N20 << ")\n"
               << ANSI_RESET;

          // -------------------------------------------------------------------------
          // Local helpers (kept minimal): merge Sam/Justin in memory, then overlay.
          // -------------------------------------------------------------------------
          auto FindXbinByEdges =
            [&](const TAxis* ax, double lo, double hi)->int
          {
            if (!ax) return -1;
            const int nb = ax->GetNbins();
            for (int ib = 1; ib <= nb; ++ib)
            {
              const double a = ax->GetBinLowEdge(ib);
              const double b = ax->GetBinUpEdge(ib);
              if (std::fabs(a - lo) < 1e-6 && std::fabs(b - hi) < 1e-6) return ib;
            }
            return -1;
          };

          auto StyleForOverlay =
            [&](TH1* h, int col)->void
          {
            if (!h) return;
            EnsureSumw2(h);
            h->SetTitle("");
            h->SetLineWidth(2);
            h->SetLineColor(col);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.00);
            h->SetMarkerColor(col);
            h->SetFillStyle(0);
            h->GetXaxis()->SetTitle("x_{J#gamma}");
            h->GetYaxis()->SetTitle("A.U.");
            h->GetXaxis()->SetRangeUser(0.0, 2.0);
            const int nb = h->GetNbinsX();
            const double integ = h->Integral(0, nb + 1);
            if (integ > 0.0) h->Scale(1.0 / integ);
          };

          auto StyleForOverlayRaw =
            [&](TH1* h, int col)->void
          {
            if (!h) return;
            EnsureSumw2(h);
            h->SetTitle("");
            h->SetLineWidth(2);
            h->SetLineColor(col);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(1.00);
            h->SetMarkerColor(col);
            h->SetFillStyle(0);
            h->GetXaxis()->SetTitle("x_{J#gamma}");
            h->GetYaxis()->SetTitle("Counts");
            h->GetXaxis()->SetRangeUser(0.0, 2.0);
          };

          auto RebinOntoTemplateAxis =
            [&](TH1*& hInOut, TH1* hTemplate)->void
          {
            if (!hInOut || !hTemplate) return;

            TH1* hReb = CloneTH1(hTemplate, TString::Format("%s_rebinnedToTemplate", hInOut->GetName()).Data());
            if (!hReb) return;

            hReb->Reset("ICES");
            hReb->SetDirectory(nullptr);
            if (hReb->GetSumw2N() == 0) hReb->Sumw2();

            const int nbIn = hInOut->GetNbinsX();
            for (int ib = 1; ib <= nbIn; ++ib)
            {
              const double x = hInOut->GetXaxis()->GetBinCenter(ib);
              const int ob = hReb->FindBin(x);
              if (ob < 1 || ob > hReb->GetNbinsX()) continue;

              const double c = hInOut->GetBinContent(ib);
              const double e = hInOut->GetBinError(ib);

              hReb->SetBinContent(ob, hReb->GetBinContent(ob) + c);

              const double oe = hReb->GetBinError(ob);
              hReb->SetBinError(ob, std::sqrt(oe*oe + e*e));
            }

            delete hInOut;
            hInOut = hReb;
          };

          auto BuildJustinMergedTH3 =
            [&](const std::string& rKey, const std::string& hPrefix, const std::string& tag)->TH3*
          {
            const std::string hname = hPrefix + rKey;

            TH3* h10 = dynamic_cast<TH3*>(dJ10->Get(hname.c_str()));
            TH3* h20 = dynamic_cast<TH3*>(dJ20->Get(hname.c_str()));
            if (!h10 && !h20) return nullptr;

            TH3* h = nullptr;

            if (h10)
            {
              h = CloneTH3(h10, TString::Format("hJustin_%s_%s", tag.c_str(), rKey.c_str()).Data());
              if (h)
              {
                h->SetDirectory(nullptr);
                if (h->GetSumw2N() == 0) h->Sumw2();
                h->Scale(jw10);
              }
            }

            if (h20)
            {
              if (!h)
              {
                h = CloneTH3(h20, TString::Format("hJustin_%s_%s", tag.c_str(), rKey.c_str()).Data());
                if (h)
                {
                  h->SetDirectory(nullptr);
                  if (h->GetSumw2N() == 0) h->Sumw2();
                  h->Scale(jw20);
                }
              }
              else
              {
                TH3* tmp = CloneTH3(h20, TString::Format("hJustin_%s_tmp20_%s", tag.c_str(), rKey.c_str()).Data());
                if (tmp)
                {
                  tmp->SetDirectory(nullptr);
                  if (tmp->GetSumw2N() == 0) tmp->Sumw2();
                  tmp->Scale(jw20);
                  h->Add(tmp);
                  delete tmp;
                }
              }
            }

            return h;
          };

          auto GetSamHistMerged =
            [&](int iPt, int jR, int kIso, int lABCD, const std::string& newName)->TH1*
          {
            const std::string nm = TString::Format("hratio_%d_%d_%d_%d_%d", iPt, jR, kIso, 0, lABCD).Data();
            TH1* h10 = dynamic_cast<TH1*>(fSam10->Get(nm.c_str()));
            TH1* h20 = dynamic_cast<TH1*>(fSam20->Get(nm.c_str()));
            if (!h10 && !h20) return nullptr;

            TH1* out = nullptr;

            if (h10)
            {
              out = CloneTH1(h10, newName);
              if (out)
              {
                out->SetDirectory(nullptr);
                if (out->GetSumw2N() == 0) out->Sumw2();
                out->Scale(sw10);
              }
            }

            if (h20)
            {
              if (!out)
              {
                out = CloneTH1(h20, newName);
                if (out)
                {
                  out->SetDirectory(nullptr);
                  if (out->GetSumw2N() == 0) out->Sumw2();
                  out->Scale(sw20);
                }
              }
              else
              {
                TH1* tmp = CloneTH1(h20, newName + "_tmp20");
                if (tmp)
                {
                  tmp->SetDirectory(nullptr);
                  if (tmp->GetSumw2N() == 0) tmp->Sumw2();
                  tmp->Scale(sw20);
                  out->Add(tmp);
                  delete tmp;
                }
              }
            }

            return out;
          };

          // Sam pT binning edges: [10,11,12,13,15,19,30]
          const std::vector<double> samPtEdges = {10,11,12,13,15,19,30};

          const int kIsoSamJES  = 1; // Sam JES-calibrated curve
          const int lABCD = 0;

          auto FindSamExactEdgeIndex =
            [&](double x)->int
          {
            for (int i = 0; i < (int)samPtEdges.size(); ++i)
            {
              if (std::fabs(samPtEdges[i] - x) < 1e-9) return i;
            }
            return -1;
          };

          auto FindSamCoveringBinIndex =
            [&](double lo, double hi)->int
          {
            for (int i = 0; i + 1 < (int)samPtEdges.size(); ++i)
            {
              const double a = samPtEdges[i];
              const double b = samPtEdges[i+1];
              if (a <= lo + 1e-9 && hi <= b + 1e-9) return i;
            }
            return -1;
          };

          auto GetSamHistForJustinBin =
            [&](double ptLo, double ptHi, int jR, int kIso_, int lABCD_, const std::string& newName, std::string& outSamLabel)->TH1*
          {
            outSamLabel = "UNKNOWN";

            const int a = FindSamExactEdgeIndex(ptLo);
            const int b = FindSamExactEdgeIndex(ptHi);

            // Exact edge match: SUM Sam bins (integration)
            if (a >= 0 && b >= 0 && b > a)
            {
              TH1* sum = nullptr;
              for (int iPt = a; iPt < b; ++iPt)
              {
                TH1* h = GetSamHistMerged(
                  iPt, jR, kIso_, lABCD_,
                  TString::Format("%s_sam_i%d", newName.c_str(), iPt).Data()
                );
                if (!h) { if (sum) delete sum; return nullptr; }

                if (!sum)
                {
                  sum = CloneTH1(h, newName);
                  if (sum)
                  {
                    sum->Reset("ICES");
                    sum->SetDirectory(nullptr);
                  }
                }
                if (sum) sum->Add(h);
                delete h;
              }

              outSamLabel = TString::Format("%.0f-%.0f", ptLo, ptHi).Data();
              return sum;
            }

            // Otherwise: use single Sam bin that COVERS [ptLo,ptHi]
            const int iCover = FindSamCoveringBinIndex(ptLo, ptHi);
            if (iCover < 0)
            {
              return nullptr;
            }

            const double sLo = samPtEdges[iCover];
            const double sHi = samPtEdges[iCover+1];
            outSamLabel = TString::Format("%.0f-%.0f", sLo, sHi).Data();

            return GetSamHistMerged(iCover, jR, kIso_, lABCD_, newName);
          };

          const auto& jes3Edges = Binning().jes3_photon_pt_bins;

          struct RMap { std::string rKey; int jIdx; double R; };
          const std::vector<RMap> rMap =
          {
            {"r04", 1, 0.4},
            {"r06", 2, 0.6}
          };

          const double jetMinPtGeV = cfgJ.jetMinPt;
          const std::string bbLabel = cfgJ.bbLabel;

          cout << ANSI_BOLD_CYN
               << "[PLOTS] Writing ONLY into: " << dirOver << "\n"
               << ANSI_RESET;

          for (const auto& rm : rMap)
          {
            const std::string& rKey = rm.rKey;
            const double R = rm.R;

            cout << ANSI_BOLD_CYN
                 << "\n[SAM VS JUSTIN] Radius: " << rKey << "  (R=" << R << ")\n"
                 << ANSI_RESET;

            TH3* hJustin3_reco = BuildJustinMergedTH3(rKey, "h_JES3_pT_xJ_alpha_", "RECO");
            if (!hJustin3_reco)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Missing Justin TH3 for " << rKey << " (h_JES3_pT_xJ_alpha_" << rKey << "). Skipping radius.\n"
                   << ANSI_RESET;
              continue;
            }

            TH3* hJustin3_tag = BuildJustinMergedTH3(rKey, "h_JES3RecoTruthTagged_pT_xJ_alpha_", "RECO_truthTaggedPhoJet");
            if (!hJustin3_tag)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Missing Justin truth-tagged TH3 for " << rKey << " (h_JES3RecoTruthTagged_pT_xJ_alpha_" << rKey << "). Will draw Sam vs Justin only.\n"
                   << ANSI_RESET;
            }

            for (int ip = 0; ip + 1 < (int)jes3Edges.size(); ++ip)
            {
              const double ptLo = jes3Edges[(std::size_t)ip];
              const double ptHi = jes3Edges[(std::size_t)ip + 1];

              const int xbin = FindXbinByEdges(hJustin3_reco->GetXaxis(), ptLo, ptHi);
              if (xbin < 1)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] Justin TH3 missing pT bin " << ptLo << "-" << ptHi << " for " << rKey << ". Skipping bin.\n"
                     << ANSI_RESET;
                continue;
              }

              TH1* hJustin = ProjectY_AtXbin_TH3(
                hJustin3_reco, xbin,
                TString::Format("hJustin_xJ_%.0f_%.0f_%s", ptLo, ptHi, rKey.c_str()).Data()
              );

              TH1* hJustinTag = nullptr;
              if (hJustin3_tag)
              {
                const int xbinTag = FindXbinByEdges(hJustin3_tag->GetXaxis(), ptLo, ptHi);
                if (xbinTag >= 1)
                {
                  hJustinTag = ProjectY_AtXbin_TH3(
                    hJustin3_tag, xbinTag,
                    TString::Format("hJustin_truthTaggedPhoJet_xJ_%.0f_%.0f_%s", ptLo, ptHi, rKey.c_str()).Data()
                  );
                }
              }

              std::string samPtLabel;
              TH1* hSam = GetSamHistForJustinBin(
                ptLo, ptHi, rm.jIdx, kIsoSamJES, lABCD,
                TString::Format("hSam_xJ_%.0f_%.0f_%s", ptLo, ptHi, rKey.c_str()).Data(),
                samPtLabel
              );

              if (!hJustin || !hSam)
              {
                if (hSam) delete hSam;
                if (hJustin) delete hJustin;
                if (hJustinTag) delete hJustinTag;
                continue;
              }

              // Sam hratio_* is typically finer-binned than Justin's xJ projection:
              // Rebin Sam onto Justin x-axis BEFORE normalization.
              RebinOntoTemplateAxis(hSam, hJustin);

                // -------------------------------------------------------------------------
                // Write TWO versions:
                //   (1) Unnormalized (Counts)  -> dirOver_raw
                //   (2) Normalized (A.U.)      -> dirOver_norm
                // Both with x-axis cut at 2.0 and a smaller legend pushed to top-right.
                // -------------------------------------------------------------------------

                auto MaxBinContent = [&](TH1* h)->double
                {
                  if (!h) return 0.0;
                  double m = 0.0;
                  const int nb = h->GetNbinsX();
                  for (int ib = 1; ib <= nb; ++ib)
                  {
                    const double v = h->GetBinContent(ib);
                    if (v > m) m = v;
                  }
                  return m;
                };

                // -------------------------
                // (1) UNNORMALIZED
                // -------------------------
                {
                  TH1* hSamRaw    = CloneTH1(hSam,    TString::Format("%s_raw", hSam->GetName()).Data());
                  TH1* hJustinRaw = CloneTH1(hJustin, TString::Format("%s_raw", hJustin->GetName()).Data());
                  TH1* hTagRaw    = (hJustinTag ? CloneTH1(hJustinTag, TString::Format("%s_raw", hJustinTag->GetName()).Data()) : nullptr);

                  if (hSamRaw)    { hSamRaw->SetDirectory(nullptr); }
                  if (hJustinRaw) { hJustinRaw->SetDirectory(nullptr); }
                  if (hTagRaw)    { hTagRaw->SetDirectory(nullptr); }

                  if (hSamRaw && hJustinRaw)
                  {
                    StyleForOverlayRaw(hSamRaw, 2);
                    StyleForOverlayRaw(hJustinRaw, 4);
                    if (hTagRaw)
                    {
                      StyleForOverlayRaw(hTagRaw, 6);
                      hTagRaw->SetMarkerStyle(24);
                    }

                    double ymaxRaw = 0.0;
                    ymaxRaw = std::max(ymaxRaw, MaxBinContent(hSamRaw));
                    ymaxRaw = std::max(ymaxRaw, MaxBinContent(hJustinRaw));
                    if (hTagRaw) ymaxRaw = std::max(ymaxRaw, MaxBinContent(hTagRaw));
                    hSamRaw->SetMaximum(ymaxRaw * 1.10);

                    TCanvas cRaw("c_SamVsJustin_raw","c_SamVsJustin_raw",900,700);
                    ApplyCanvasMargins1D(cRaw);

                    hSamRaw->Draw("E1");
                    hJustinRaw->Draw("E1 same");
                    if (hTagRaw) hTagRaw->Draw("E1 same");

                    // Optional cut lines ONLY for the 13-15 GeV bin (as before)
                    TLine* lAbs  = nullptr;
                    TLine* lFull = nullptr;

                    const bool is1315 = (std::fabs(ptLo - 13.0) < 1e-6 && std::fabs(ptHi - 15.0) < 1e-6);
                    if (is1315)
                    {
                      const double xAbs  = jetMinPtGeV / ptHi;
                      const double xFull = jetMinPtGeV / ptLo;

                      lAbs  = new TLine(xAbs,  0.0, xAbs,  ymaxRaw * 1.10);
                      lFull = new TLine(xFull, 0.0, xFull, ymaxRaw * 1.10);

                      lAbs->SetLineColor(4);
                      lAbs->SetLineWidth(2);
                      lAbs->SetLineStyle(2);

                      lFull->SetLineColor(2);
                      lFull->SetLineWidth(2);
                      lFull->SetLineStyle(2);

                      lAbs->Draw("same");
                      lFull->Draw("same");
                    }

                    // Smaller legend, pushed to top-right
                    TLegend leg(0.62, 0.82, 0.90, 0.95);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.022);
                    leg.SetFillStyle(0);
                    leg.SetBorderSize(0);
                    leg.AddEntry(hSamRaw,    TString::Format("Sam's RECO (R = %.1f)", R).Data(), "ep");
                    leg.AddEntry(hJustinRaw, TString::Format("Justin's RECO w/ matched cuts (R = %.1f)", R).Data(), "ep");
                    if (hTagRaw)
                    {
                      leg.AddEntry(hTagRaw, TString::Format("Justin's RECO (#gamma^{truth} + jet^{truth} tag) (R = %.1f)", R).Data(), "ep");
                    }
                    leg.Draw();

                    // Note: aligned under the legend block
                    TLatex tNote;
                    tNote.SetNDC(true);
                    tNote.SetTextFont(42);
                    tNote.SetTextAlign(13);
                    tNote.SetTextSize(0.032);
                    tNote.DrawLatex(0.62, 0.79, "Photon 10 + 20 GeV #gamma+Jet MC");

                    // Info block: middle RHS
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.032);
                    t.SetTextAlign(12);
                    t.DrawLatex(0.60, 0.62, TString::Format("p_{T}^{#gamma}: %.0f-%.0f GeV", ptLo, ptHi).Data());
                    t.DrawLatex(0.60, 0.57, TString::Format("Sam p_{T}^{#gamma} used: %s GeV", samPtLabel.c_str()).Data());
                    t.DrawLatex(0.60, 0.52, TString::Format("p_{T}^{jet,min} = %.0f GeV", jetMinPtGeV).Data());
                    t.DrawLatex(0.60, 0.47, TString::Format("Back-to-back: %s", bbLabel.c_str()).Data());

                    TLegend* legCuts = nullptr;
                    if (is1315 && lAbs && lFull)
                    {
                      legCuts = new TLegend(0.60, 0.32, 0.88, 0.43);
                      legCuts->SetTextFont(42);
                      legCuts->SetTextSize(0.030);
                      legCuts->SetFillStyle(0);
                      legCuts->SetBorderSize(0);
                      legCuts->AddEntry(lAbs,  "x_{J,min}^{abs} = p_{T}^{jet,min}/p_{T,max}^{#gamma}", "l");
                      legCuts->AddEntry(lFull, "x_{J,min}^{full} = p_{T}^{jet,min}/p_{T,min}^{#gamma}", "l");
                      legCuts->Draw();
                    }

                    const std::string outNameRaw =
                      TString::Format("overlay_SamVsJustin_JES3_RECO_RAW_pTgamma_%.0f_%.0f_Sam_%s_%s.png",
                        ptLo, ptHi, samPtLabel.c_str(), rKey.c_str()).Data();

                    cRaw.Modified();
                    cRaw.Update();
                    SaveCanvas(cRaw, JoinPath(dirOver_raw, outNameRaw));

                    if (legCuts) delete legCuts;
                    if (lAbs)  delete lAbs;
                    if (lFull) delete lFull;
                  }

                  if (hSamRaw) delete hSamRaw;
                  if (hJustinRaw) delete hJustinRaw;
                  if (hTagRaw) delete hTagRaw;
                }

                // -------------------------
                // (2) NORMALIZED
                // -------------------------
                {
                  StyleForOverlay(hSam, 2);
                  StyleForOverlay(hJustin, 4);
                  if (hJustinTag)
                  {
                    StyleForOverlay(hJustinTag, 6);
                    hJustinTag->SetMarkerStyle(24);
                  }

                  double ymax = 0.0;
                  ymax = std::max(ymax, MaxBinContent(hSam));
                  ymax = std::max(ymax, MaxBinContent(hJustin));
                  if (hJustinTag) ymax = std::max(ymax, MaxBinContent(hJustinTag));
                  hSam->SetMaximum(ymax * 1.10);

                  TCanvas c("c_SamVsJustin","c_SamVsJustin",900,700);
                  ApplyCanvasMargins1D(c);

                  hSam->Draw("E1");
                  hJustin->Draw("E1 same");
                  if (hJustinTag) hJustinTag->Draw("E1 same");

                  // Optional cut lines ONLY for the 13-15 GeV bin (as before)
                  TLine* lAbs  = nullptr;
                  TLine* lFull = nullptr;

                  const bool is1315 = (std::fabs(ptLo - 13.0) < 1e-6 && std::fabs(ptHi - 15.0) < 1e-6);
                  if (is1315)
                  {
                    const double xAbs  = jetMinPtGeV / ptHi;
                    const double xFull = jetMinPtGeV / ptLo;

                    lAbs  = new TLine(xAbs,  0.0, xAbs,  ymax * 1.10);
                    lFull = new TLine(xFull, 0.0, xFull, ymax * 1.10);

                    lAbs->SetLineColor(4);
                    lAbs->SetLineWidth(2);
                    lAbs->SetLineStyle(2);

                    lFull->SetLineColor(2);
                    lFull->SetLineWidth(2);
                    lFull->SetLineStyle(2);

                    lAbs->Draw("same");
                    lFull->Draw("same");
                  }

                  // Smaller legend, pushed to top-right
                  TLegend leg(0.62, 0.82, 0.90, 0.95);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.022);
                  leg.SetFillStyle(0);
                  leg.SetBorderSize(0);
                  leg.AddEntry(hSam,    TString::Format("Sam's RECO (R = %.1f)", R).Data(), "ep");
                  leg.AddEntry(hJustin, TString::Format("Justin's RECO w/ matched cuts (R = %.1f)", R).Data(), "ep");
                  if (hJustinTag)
                  {
                    leg.AddEntry(hJustinTag, TString::Format("Justin's RECO (#gamma^{truth} + jet^{truth} tag) (R = %.1f)", R).Data(), "ep");
                  }
                  leg.Draw();

                  // Note: aligned under the legend block
                  TLatex tNote;
                  tNote.SetNDC(true);
                  tNote.SetTextFont(42);
                  tNote.SetTextAlign(13);
                  tNote.SetTextSize(0.032);
                  tNote.DrawLatex(0.62, 0.79, "Photon 10 + 20 GeV #gamma+Jet MC");

                  // Info block: middle RHS
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.032);
                  t.SetTextAlign(12);
                  t.DrawLatex(0.60, 0.62, TString::Format("p_{T}^{#gamma}: %.0f-%.0f GeV", ptLo, ptHi).Data());
                  t.DrawLatex(0.60, 0.57, TString::Format("Sam p_{T}^{#gamma} used: %s GeV", samPtLabel.c_str()).Data());
                  t.DrawLatex(0.60, 0.52, TString::Format("p_{T}^{jet,min} = %.0f GeV", jetMinPtGeV).Data());
                  t.DrawLatex(0.60, 0.47, TString::Format("Back-to-back: %s", bbLabel.c_str()).Data());

                  TLegend* legCuts = nullptr;
                  if (is1315 && lAbs && lFull)
                  {
                    legCuts = new TLegend(0.60, 0.32, 0.88, 0.43);
                    legCuts->SetTextFont(42);
                    legCuts->SetTextSize(0.030);
                    legCuts->SetFillStyle(0);
                    legCuts->SetBorderSize(0);
                    legCuts->AddEntry(lAbs,  "x_{J,min}^{abs} = p_{T}^{jet,min}/p_{T,max}^{#gamma}", "l");
                    legCuts->AddEntry(lFull, "x_{J,min}^{full} = p_{T}^{jet,min}/p_{T,min}^{#gamma}", "l");
                    legCuts->Draw();
                  }

                  const std::string outName =
                    TString::Format("overlay_SamVsJustin_JES3_RECO_NORM_pTgamma_%.0f_%.0f_Sam_%s_%s.png",
                      ptLo, ptHi, samPtLabel.c_str(), rKey.c_str()).Data();

                  c.Modified();
                  c.Update();
                  SaveCanvas(c, JoinPath(dirOver_norm, outName));

                  if (legCuts) delete legCuts;
                  if (lAbs)  delete lAbs;
                  if (lFull) delete lFull;
                }

                delete hSam;
                delete hJustin;
                if (hJustinTag) delete hJustinTag;
            }

            delete hJustin3_reco;
            if (hJustin3_tag) delete hJustin3_tag;
          }

          fJ10->Close();
          fJ20->Close();

          fSam10->Close();
          fSam20->Close();

          cout << ANSI_BOLD_GRN
               << "\n[OK] Sam-vs-Justin overlays complete. All PNGs written under:\n"
               << "     " << dirOver << "\n"
               << ANSI_RESET;
      }

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

                            const auto& cfgDef = DefaultSim10and20Config();
                            const double jetPtMin_GeV = cfgDef.jetMinPt;
                            const string bbLabel = cfgDef.bbLabel;

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

                            TLegend* leg = new TLegend(0.7, 0.28, 0.95, 0.45);
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
                    const int nRows = 2;
                    const int perPage = nCols * nRows;

                    const int startBinForTable = 1;
                    const int nTableBins = std::min(perPage, nPt);

                    TCanvas canTbl(
                      TString::Format("c_tbl_%s_dataVsSim", rKey.c_str()).Data(),
                      "c_tbl_dataVsSim", 1500, 900
                    );
                    canTbl.Divide(nCols, nRows, 0.001, 0.001);

                    std::vector<TH1*> keep;
                    keep.reserve(2 * nTableBins);

                    for (int k = 0; k < nTableBins; ++k)
                    {
                      const int ib = startBinForTable + k;
                      canTbl.cd(k + 1);

                      gPad->SetLeftMargin(0.14);
                      gPad->SetRightMargin(0.05);
                      gPad->SetTopMargin(0.12);
                      gPad->SetBottomMargin(0.14);

                      const double ptMinGamma = H.hReco_xJ->GetXaxis()->GetBinLowEdge(ib);
                      const double ptMaxGamma = H.hReco_xJ->GetXaxis()->GetBinUpEdge(ib);

                      TH1* hDatRaw = ProjectY_AtXbin_AndAlphaMax_TH3(
                        H.hReco_xJ, ib, H.hReco_xJ->GetZaxis()->GetXmax(),
                        TString::Format("h_tbl_dat_%s_%d", rKey.c_str(), ib).Data()
                      );
                      TH1* hSimRaw = ProjectY_AtXbin_AndAlphaMax_TH3(
                        hSim3, ib, hSim3->GetZaxis()->GetXmax(),
                        TString::Format("h_tbl_sim_%s_%d", rKey.c_str(), ib).Data()
                      );

                      if (!hDatRaw || !hSimRaw) { if (hDatRaw) delete hDatRaw; if (hSimRaw) delete hSimRaw; continue; }

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

                        const auto& cfgDef = DefaultSim10and20Config();
                        const double jetPtMin_GeV = cfgDef.jetMinPt;
                        const string bbLabel = cfgDef.bbLabel;

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

                        TLine* lnAbsDrawn  = nullptr;
                        TLine* lnFullDrawn = nullptr;

                        if (xAbs > 0.0)  lnAbsDrawn  = (TLine*) lnAbs->DrawClone("same");
                        if (xFull > 0.0) lnFullDrawn = (TLine*) lnFull->DrawClone("same");

                        delete lnFull;
                        delete lnAbs;

                        TLegend* leg = new TLegend(0.52, 0.70, 0.95, 0.90);
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.032);

                        leg->AddEntry(hDatRaw, "DATA (reco)", "ep");
                        leg->AddEntry(hSimRaw, "SIM (reco)",  "ep");
                        if (lnAbsDrawn)
                          leg->AddEntry(lnAbsDrawn,
                            TString::Format("x_{J, min}^{abs} = #frac{%.0f}{p_{T, max}^{#gamma}} = %.3f", jetPtMin_GeV, xAbs),
                            "l");
                        if (lnFullDrawn)
                          leg->AddEntry(lnFullDrawn,
                            TString::Format("x_{J, min}^{full} = #frac{%.0f}{p_{T, min}^{#gamma}} = %.3f", jetPtMin_GeV, xFull),
                            "l");

                        leg->DrawClone();
                        delete leg;

                        {
                          TLatex tCuts;
                          tCuts.SetNDC(true);
                          tCuts.SetTextFont(42);
                          tCuts.SetTextAlign(33);
                          tCuts.SetTextSize(0.038);
                          tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data());
                          tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                        }

                        TLatex ttl;
                        ttl.SetNDC(true);
                        ttl.SetTextFont(42);
                        ttl.SetTextSize(0.052);
                        ttl.DrawLatex(0.12, 0.94,
                          TString::Format("RECO x_{J#gamma} (DATA vs SIM), p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                            ptMinGamma, ptMaxGamma, R).Data());

                        keep.push_back(hDatRaw);
                        keep.push_back(hSimRaw);
                    }

                    SaveCanvas(canTbl, JoinPath(dirOv, "table3x2_overlay_integratedAlpha_overlayedWithSim.png"));

                    for (auto* h1 : keep) delete h1;
                  }
                }

                // -------------------------------------------------------------------------
                // 3x3 tables + text summary (existing behavior preserved)
                // -------------------------------------------------------------------------

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
                      "c_tbl_xJ", 1500, (perPage == 6 ? 900 : 1200)
                    );

                      c.Divide(nCols, nRows, 0.001, 0.001);

                      std::vector<TObject*> keep;
                      keep.reserve(perPage * 4);

                      for (int k = 0; k < perPage; ++k)
                  {
                    const int ib = start + k;
                    c.cd(k+1);

                    gPad->SetLeftMargin(0.14);
                    gPad->SetRightMargin(0.05);
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

                      const auto& cfgDef = DefaultSim10and20Config();
                      const double jetPtMin_GeV = cfgDef.jetMinPt;

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

                      TLegend* leg = new TLegend(0.44, 0.67, 0.96, 0.90);
                      leg->SetBorderSize(0);
                      leg->SetFillStyle(0);
                      leg->SetTextFont(42);
                      leg->SetTextSize(0.038);

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
                        tCuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", cfgDef.bbLabel.c_str()).Data());
                        tCuts.DrawLatex(0.92, 0.54, TString::Format("p_{T}^{jet} > %.0f GeV", jetPtMin_GeV).Data());
                      }

                      TLatex ttl;
                      ttl.SetNDC(true);
                      ttl.SetTextFont(42);
                      ttl.SetTextSize(0.052);
                      ttl.DrawLatex(0.14, 0.94,
                        TString::Format("RECO x_{J#gamma}, p_{T}^{#gamma} = %.0f - %.0f GeV, R = %.1f",
                          ptMin, ptMax, R).Data());

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
              Make3x3Table_xJ_FromTH3(H.hReco_xJ,               D.dirXJProjReco,                   "RECO",  false);
              Make3x3Table_xJ_FromTH3(H.hRecoTruthPhoTagged_xJ, D.dirXJProjRecoTruthPhoTagged,     "RECO",  false);
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
                // NEW: Overlays (shape), integrated over alpha
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
                  if (!hBlack || !hRed) return;

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

                      // Legend: shift slightly left and make only slightly smaller
                      TLegend leg(0.55, 0.75, 0.85, 0.90);
                      leg.SetTextFont(42);
                      leg.SetTextSize(0.038);
                      leg.SetFillStyle(0);
                      leg.SetBorderSize(0);
                      leg.AddEntry(hA, legBlack.c_str(), "ep");
                      leg.AddEntry(hB, legRed.c_str(),   "ep");
                      leg.Draw();

                      // Clean, minimal top-left header (no "Overlay (shape)", no jet-match line, no "Integrated over α")
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
                      const std::string bbLabel = DefaultSim10and20Config().bbLabel;
                      const double jetMinPt = DefaultSim10and20Config().jetMinPt;

                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextAlign(13); // left-top
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

                    SaveCanvas(c, outPng);

                    delete hA;
                    delete hB;
                  }

                    // --- 2x3 table page of overlays (shape): last 6 pT bins (skip first pT bin) ---
                    const int nCols   = 3;
                    const int nRows   = 2;
                    const int perPage = nCols * nRows;  // 6

                    const int firstBinToTable = 2;  // skip the first pT bin
                    const int startBin = (ovTag == "RECO_vs_RECO_truthTaggedPhoJet")
                      ? 1
                      : std::max(firstBinToTable, nPtOv - perPage + 1);

                    if (startBin <= nPtOv)
                    {
                      const int page = 1;

                      TCanvas c(
                        TString::Format("c_tbl_ov_%s_%s_p%d", ovTag.c_str(), rKey.c_str(), page).Data(),
                        "c_tbl_ov", 1500, 900
                      );
                      c.Divide(nCols, nRows, 0.001, 0.001);

                      vector<TH1*> keep;
                      keep.reserve(2 * perPage);

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

                          const std::string bbLabel = DefaultSim10and20Config().bbLabel;
                          const double jetMinPt     = DefaultSim10and20Config().jetMinPt;
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

                      const string outName = "table2x3_overlay_shape.png";

                      SaveCanvas(c, JoinPath(dirOvBase, outName));

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
                  if (!hReco || !hTruth || !hRecoTruth) return;

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

                    // Kinematic floor line (requested)
                    DrawKinematicFloorLine(xFloor, 0.0, yMaxPlot);

                    // Small in-pad label for the floor
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

                    DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);

                    vector<string> lines = headerLines;
                    lines.push_back(TString::Format("p_{T}^{#gamma}: %s  (R=%.1f)", ptLab.c_str(), R).Data());
                    lines.push_back("Integrated over #alpha");
                    DrawLatexLines(0.14, 0.84, lines, 0.030, 0.040);

                    SaveCanvas(c, outPng);

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

                              const auto& simCfgMap = Sim10and20Configs();
                              const auto  itSimCfg  = simCfgMap.find(DefaultSimSampleKey());
                              const Sim10and20Config* simCfg = (itSimCfg != simCfgMap.end() ? &itSimCfg->second : nullptr);

                              DrawLatexLines(0.14, 0.92,
                                {TString::Format("Lead recoil-jet match efficiency (R = %.1f)", Rval).Data()},
                                0.030, 0.040
                              );
                              DrawLatexLines(0.14, 0.86, DefaultHeaderLines(ds), 0.030, 0.040);

                              if (simCfg)
                              {
                                TLatex tCuts;
                                tCuts.SetNDC(true);
                                tCuts.SetTextFont(42);
                                tCuts.SetTextAlign(33);
                                tCuts.SetTextSize(0.028);
                                tCuts.DrawLatex(0.92, 0.78, TString::Format("|#Delta#phi(#gamma,jet)| > %s", simCfg->bbLabel.c_str()).Data());
                                tCuts.DrawLatex(0.92, 0.70, TString::Format("Reco (p_{T}^{jet} > %.0f GeV)", simCfg->jetMinPt).Data());
                              }


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

                                  // Default SIM cuts (from kDefaultSimSampleKey) — print UNDER the title block (2 lines)
                                  const auto& cfgDef = DefaultSim10and20Config();
                                  DrawLatexLines(0.14, 0.80,
                                    {
                                      TString::Format("#Delta#phi(#gamma,jet) > %s", cfgDef.bbLabel.c_str()).Data(),
                                      TString::Format("p_{T}^{jet} > %.0f GeV", cfgDef.jetMinPt).Data()
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


        // =============================================================================
        // Section 5G: (pTgamma,eta,phi) maps: photon, recoilJet1, and balance profile
        // =============================================================================
        void RunMapsQA(Dataset& ds)
        {
          cout << ANSI_BOLD_CYN << "\n==============================\n"
               << "[SECTION 5G] Maps (pTgamma,#eta,#phi) (" << ds.label << ")\n"
               << "==============================" << ANSI_RESET << "\n";

          string outDir = ds.isSim
              ? JoinPath(ds.outBase, "RecoilJetQA/Maps")
              : JoinPath(ds.outBase, "baselineData/RecoilJetQA/Maps");
          EnsureDir(outDir);

          // Photon maps (tight+isolated photons) if available
          if (TH3* hPho = GetObj<TH3>(ds, "h_Pho3_pT_eta_phi_tightIso", true, true, true))
          {
            const string phoDir = JoinPath(outDir, "Photon");
            EnsureDir(phoDir);

            // integrated eta-phi
            {
              TH3* hc = CloneTH3(hPho, "pho3_clone");
              hc->GetXaxis()->SetRange(1, hc->GetXaxis()->GetNbins());
              TH2* h2 = dynamic_cast<TH2*>(hc->Project3D("zy"));

              if (h2)
              {
                h2->SetDirectory(nullptr);
                vector<string> lines = {"Photon #eta-#phi map (tightIso)", "Integrated over p_{T}^{#gamma}"};
                DrawAndSaveTH2_Common(ds, h2,
                  JoinPath(phoDir, "pho_etaPhi_integrated.png"),
                  "#eta^{#gamma}", "#phi^{#gamma}", "Counts", lines, false);
                delete h2;
              }
              delete hc;
            }

            // per pT bin projections (use axis binning, not PtBins hardcoding)
            const int nPt = hPho->GetXaxis()->GetNbins();
            for (int ib = 1; ib <= nPt; ++ib)
            {
              const string ptLab = AxisBinLabel(hPho->GetXaxis(), ib, "GeV", 0);
              TH2* h2 = ProjectYZ_AtXbin_TH3(hPho, ib, TString::Format("pho_etaPhi_%d", ib).Data());
              if (!h2) continue;

              vector<string> lines = {"Photon #eta-#phi map (tightIso)", TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()};
              DrawAndSaveTH2_Common(ds, h2,
                JoinPath(phoDir, TString::Format("pho_etaPhi_pTbin%d.png", ib).Data()),
                "#eta^{#gamma}", "#phi^{#gamma}", "Counts", lines, false);
              delete h2;
            }
          }

          // RecoilJet1 maps + Balance profiles per rKey (SIM and/or DATA)
          for (const auto& rKey : kRKeys)
          {
            const double R = RFromKey(rKey);
            const double etaFidAbs = FidEtaAbsFromKey(rKey);

            const string rDir = JoinPath(outDir, rKey);
            const string jetDir = JoinPath(rDir, "RecoilJet1");
            const string balDir = JoinPath(rDir, "Balance");
            EnsureDir(rDir);
            EnsureDir(jetDir);
            EnsureDir(balDir);

            // TH3: h_Jet13_pTgamma_eta_phi_recoilJet1_rXX
            if (TH3* hJ = GetObj<TH3>(ds, "h_Jet13_pTgamma_eta_phi_recoilJet1_" + rKey, true, true, true))
            {
              // integrated eta-phi
              {
                TH3* hc = CloneTH3(hJ, "jet13_clone");
                hc->GetXaxis()->SetRange(1, hc->GetXaxis()->GetNbins());
                  TH2* h2 = dynamic_cast<TH2*>(hc->Project3D("zy"));
                  if (h2)
                  {
                    h2->SetDirectory(nullptr);
                    vector<string> lines = {
                      "RecoilJet1 #eta-#phi map",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      "Integrated over p_{T}^{#gamma}"
                    };
                    DrawAndSaveTH2_Common(ds, h2,
                      JoinPath(jetDir, "recoilJet1_etaPhi_integrated.png"),
                      "#eta^{jet1}", "#phi^{jet1}", "Counts", lines,
                      false, true, etaFidAbs);
                    delete h2;
                  }

                delete hc;
              }

              // per pT bin
              const int nPt = hJ->GetXaxis()->GetNbins();
              for (int ib = 1; ib <= nPt; ++ib)
              {
                const string ptLab = AxisBinLabel(hJ->GetXaxis(), ib, "GeV", 0);
                TH2* h2 = ProjectYZ_AtXbin_TH3(hJ, ib, TString::Format("jet13_etaPhi_%s_%d", rKey.c_str(), ib).Data());
                if (!h2) continue;

                vector<string> lines = {
                  "RecoilJet1 #eta-#phi map",
                  TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                  TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
                };

                DrawAndSaveTH2_Common(ds, h2,
                  JoinPath(jetDir, TString::Format("recoilJet1_etaPhi_pTbin%d.png", ib).Data()),
                  "#eta^{jet1}", "#phi^{jet1}", "Counts", lines,
                  false, true, etaFidAbs);

                delete h2;
              }
            }

            // TProfile3D: p_Balance3_pTgamma_eta_phi_rXX (mean xJ per (pT,eta,phi) bin)
            if (TProfile3D* pB = GetObj<TProfile3D>(ds, "p_Balance3_pTgamma_eta_phi_" + rKey, true, true, true))
            {
              // integrated over pT: take full x range, project to eta-phi profile
              {
                TProfile3D* pc = (TProfile3D*)pB->Clone("bal3_clone");
                pc->SetDirectory(nullptr);
                pc->GetXaxis()->SetRange(1, pc->GetXaxis()->GetNbins());
                  TProfile2D* p2 = dynamic_cast<TProfile2D*>(pc->Project3DProfile("zy"));
                  if (p2)
                  {
                    p2->SetDirectory(nullptr);
                    vector<string> lines = {
                      "Balance map: <x_{J}> in (#eta,#phi)",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      "Integrated over p_{T}^{#gamma}"
                    };
                    DrawAndSaveTH2_Common(ds, (TH2*)p2,
                      JoinPath(balDir, "balance_meanxJ_etaPhi_integrated.png"),
                      "#eta^{jet1}", "#phi^{jet1}", "<x_{J}>", lines,
                      false, true, etaFidAbs);
                    delete p2;
                  }

                delete pc;
              }

              // per pT bin: profile2D
              const int nPt = pB->GetXaxis()->GetNbins();
              for (int ib = 1; ib <= nPt; ++ib)
              {
                const string ptLab = AxisBinLabel(pB->GetXaxis(), ib, "GeV", 0);
                TProfile2D* p2 = ProjectYZ_AtXbin_Profile3D(pB, ib, TString::Format("bal_etaPhi_%s_%d", rKey.c_str(), ib).Data());
                if (!p2) continue;

                vector<string> lines = {
                  "Balance map: <x_{J}> in (#eta,#phi)",
                  TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                  TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data()
                };

                DrawAndSaveTH2_Common(ds, (TH2*)p2,
                  JoinPath(balDir, TString::Format("balance_meanxJ_etaPhi_pTbin%d.png", ib).Data()),
                  "#eta^{jet1}", "#phi^{jet1}", "<x_{J}>", lines,
                  false, true, etaFidAbs);

                delete p2;
              }
            }
          }
        }

        // =============================================================================
        // Section 5H: Unfolding QA for (pTgamma, xJ) and response matrices
        // =============================================================================
        void RunUnfoldingQA(Dataset& ds)
        {
            cout << ANSI_BOLD_CYN << "\n==============================\n"
                 << "[SECTION 5H] Unfolding QA (pTgamma,xJ) (" << ds.label << ")\n"
                 << "==============================" << ANSI_RESET << "\n";

            string outDir = ds.isSim
              ? JoinPath(ds.outBase, "RecoilJetQA/Unfolding")
              : JoinPath(ds.outBase, "unfolding");

            EnsureDir(outDir);

            // --------------------------------------------------------------------------
            // Helpers (refactor-only: preserve filenames, directory structure, logic)
            // --------------------------------------------------------------------------

            struct UnfoldHists
            {
              // Core xJ unfolding hists
              TH2* hReco   = nullptr;
              TH2* hTruth  = nullptr;
              TH2* hFakes  = nullptr;
              TH2* hMisses = nullptr;
              TH2* hResp   = nullptr;

              // Δphi unfolding (inclusive)
              TH2* hRecoDphi  = nullptr;
              TH2* hTruthDphi = nullptr;
            };

            auto LoadUnfoldHists =
              [&](const string& rKey)->UnfoldHists
            {
              UnfoldHists H;

              // RECO-side inputs exist in SIM and DATA
              H.hReco    = GetObj<TH2>(ds, "h2_unfoldReco_pTgamma_xJ_incl_" + rKey, true, true, true);
              H.hRecoDphi = GetObj<TH2>(ds, "h2_unfoldReco_pTgamma_dphi_incl_" + rKey, true, true, true);

              // SIM-only inputs: do NOT try to load these in DATA mode
              if (ds.isSim)
              {
                H.hTruth   = GetObj<TH2>(ds, "h2_unfoldTruth_pTgamma_xJ_incl_" + rKey, true, true, true);
                H.hFakes   = GetObj<TH2>(ds, "h2_unfoldRecoFakes_pTgamma_xJ_incl_" + rKey, true, true, true);
                H.hMisses  = GetObj<TH2>(ds, "h2_unfoldTruthMisses_pTgamma_xJ_incl_" + rKey, true, true, true);
                H.hResp    = GetObj<TH2>(ds, "h2_unfoldResponse_pTgamma_xJ_incl_" + rKey, true, true, true);

                H.hTruthDphi = GetObj<TH2>(ds, "h2_unfoldTruth_pTgamma_dphi_incl_" + rKey, true, true, true);
              }

              return H;
            };

            auto HasAnyUnfold =
              [&](const UnfoldHists& H)->bool
            {
              return (H.hReco || H.hTruth || H.hResp || H.hRecoDphi || H.hTruthDphi);
            };

            // 2D map saver (kept identical to your lambda save2D)
            auto Save2D =
              [&](const string& outBaseDir,
                  TH2* h,
                  const string& fname,
                  const string& xTitle,
                  const string& yTitle,
                  const string& zTitle,
                  const vector<string>& extra,
                  bool logz)
            {
              if (!h) return;
              TH2* hc = CloneTH2(h, fname + "_clone");
              if (!hc) return;
              DrawAndSaveTH2_Common(ds, hc,
                JoinPath(outBaseDir, fname),
                xTitle, yTitle, zTitle,
                extra, logz);
              delete hc;
            };

            // --------------------------------------------------------------------------
            // Part A (runs for SIM + DATA): core xJ unfolding 2D inputs + response visuals
            //   - preserves all filenames you had for these products
            // --------------------------------------------------------------------------
            auto RunCoreUnfolding2D =
              [&](const string& rKey, double R, const string& rOut, const UnfoldHists& H)
            {
              // Core 2D maps (presentation-ready)
              Save2D(rOut, H.hReco, "unfold_reco_pTgamma_vs_xJ.png",
                     "p_{T}^{#gamma,reco} [GeV]", "x_{J#gamma}^{reco}", "Counts",
                     {string("Unfold input: RECO counts"), rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              Save2D(rOut, H.hTruth, "unfold_truth_pTgamma_vs_xJ.png",
                     "p_{T}^{#gamma,truth} [GeV]", "x_{J#gamma}^{truth}", "Counts",
                     {string("Unfold truth: TRUTH counts"), rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              Save2D(rOut, H.hFakes, "unfold_recoFakes_pTgamma_vs_xJ.png",
                     "p_{T}^{#gamma,reco} [GeV]", "x_{J#gamma}^{reco}", "Counts (fakes)",
                     {string("Unfold input: RECO fakes"), rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              // Response matrix (SCAT) (kept exactly as your block)
              if (H.hResp)
              {
                TH2* hc = CloneTH2(H.hResp, "resp_scatter_clone");
                if (hc)
                {
                  const int oldOptStat = gStyle->GetOptStat();
                  gStyle->SetOptStat(0);
                  hc->SetStats(0);

                  TCanvas c("c_resp_scatter","c_resp_scatter",950,780);
                  ApplyCanvasMargins2D(c);
                  c.SetLogz(1);

                  hc->SetTitle("");
                  hc->GetXaxis()->SetTitle("global bin (truth: p_{T}^{#gamma}, x_{J})");
                  hc->GetYaxis()->SetTitle("global bin (reco: p_{T}^{#gamma}, x_{J})");
                  hc->GetZaxis()->SetTitle("Counts");

                  // required for log-z (must be > 0)
                  hc->SetMinimum(0.5);

                  // Use the Z palette (this is what gives you the yellow/high-count regions)
                  hc->Draw("COLZ");

                  DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                  DrawLatexLines(0.14, 0.84,
                    { "Unfold response matrix (global-bin indexing)",
                      rKey + TString::Format(" (R=%.1f)", R).Data() },
                    0.030, 0.040);

                  SaveCanvas(c, JoinPath(rOut, "unfold_response_globalTruth_vs_globalReco_SCAT.png"));

                  gStyle->SetOptStat(oldOptStat);

                  delete hc;
                }
              }
            };

            // --------------------------------------------------------------------------
            // Part B (SIM-centric): Δphi unfolding inputs + projections (only runs if present)
            // --------------------------------------------------------------------------
            auto RunDeltaPhiUnfoldingQA =
              [&](const string& rKey, double R, const string& rOut, const UnfoldHists& H)
            {
              if (!(H.hRecoDphi || H.hTruthDphi)) return;
                
                const string dphiDir = JoinPath(rOut, "deltaPhiInclusive");
                const string dphiRecoDir  = JoinPath(dphiDir, "RECO");
                const string dphiTruthDir = JoinPath(dphiDir, "TRUTH");
                const string dphiProjDir  = JoinPath(dphiDir, "projections");
                EnsureDir(dphiDir);
                EnsureDir(dphiRecoDir);
                EnsureDir(dphiProjDir);

                // Only create TRUTH folders in SIM mode (prevents empty TRUTH dirs in DATA)
                if (ds.isSim && H.hTruthDphi) EnsureDir(dphiTruthDir);


              // 2D maps (presentation-ready) — NOTE: uses Save2D with rOut path to preserve your output name
              Save2D(rOut, H.hRecoDphi,  "unfold_reco_pTgamma_vs_absDphi.png",
                     "p_{T}^{#gamma,reco} [GeV]",  "|#Delta#phi(#gamma,jet)| [rad]", "Counts",
                     {string("Unfold input: RECO |#Delta#phi| (inclusive over recoil jets)"),
                      rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              Save2D(rOut, H.hTruthDphi, "unfold_truth_pTgamma_vs_absDphi.png",
                     "p_{T}^{#gamma,truth} [GeV]", "|#Delta#phi(#gamma,jet)| [rad]", "Counts",
                     {string("Unfold truth: TRUTH |#Delta#phi| (inclusive over truth recoil jets)"),
                      rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              // Helpers: ProjectionY + visible-bin normalization
              auto ProjY_AllX = [&](TH2* h2, const string& newName)->TH1D*
              {
                if (!h2) return nullptr;
                const int nx = h2->GetXaxis()->GetNbins();
                TH1D* p = h2->ProjectionY(newName.c_str(), 1, nx);
                if (p) p->SetDirectory(nullptr);
                return p;
              };

              auto ProjY_AtXbin = [&](TH2* h2, int xbin, const string& newName)->TH1D*
              {
                if (!h2) return nullptr;
                TH1D* p = h2->ProjectionY(newName.c_str(), xbin, xbin);
                if (p) p->SetDirectory(nullptr);
                return p;
              };

              auto NormalizeVisible = [&](TH1* h)
              {
                if (!h) return;
                const int nb = h->GetNbinsX();
                const double integral = h->Integral(1, nb);
                if (integral > 0.0) h->Scale(1.0 / integral);
              };

              // 1D integrated over pTγ
              if (H.hRecoDphi)
              {
                TH1D* p = ProjY_AllX(H.hRecoDphi, TString::Format("p_absDphi_reco_all_%s_%s", ds.label.c_str(), rKey.c_str()).Data());
                if (p)
                {
                  vector<string> lines = {
                    "RECO inclusive |#Delta#phi(#gamma,jet)| (all p_{T}^{#gamma})",
                    TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data()
                  };

                  DrawAndSaveTH1_Common(ds, p,
                    JoinPath(dphiRecoDir, "absDphi_integrated_counts.png"),
                    "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1");

                  TH1* ps = CloneTH1(p, TString::Format("p_absDphi_reco_all_shape_%s_%s", ds.label.c_str(), rKey.c_str()).Data());
                  NormalizeVisible(ps);

                  DrawAndSaveTH1_Common(ds, ps,
                    JoinPath(dphiRecoDir, "absDphi_integrated_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", lines, false, false, 0.0, "E1");

                  delete ps;
                  delete p;
                }
              }

              if (H.hTruthDphi)
              {
                TH1D* p = ProjY_AllX(H.hTruthDphi, TString::Format("p_absDphi_truth_all_%s_%s", ds.label.c_str(), rKey.c_str()).Data());
                if (p)
                {
                  vector<string> lines = {
                    "TRUTH inclusive |#Delta#phi(#gamma,jet)| (all p_{T}^{#gamma})",
                    TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data()
                  };

                  DrawAndSaveTH1_Common(ds, p,
                    JoinPath(dphiTruthDir, "absDphi_integrated_counts.png"),
                    "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1");

                  TH1* ps = CloneTH1(p, TString::Format("p_absDphi_truth_all_shape_%s_%s", ds.label.c_str(), rKey.c_str()).Data());
                  NormalizeVisible(ps);

                  DrawAndSaveTH1_Common(ds, ps,
                    JoinPath(dphiTruthDir, "absDphi_integrated_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", lines, false, false, 0.0, "E1");

                  delete ps;
                  delete p;
                }
              }

              // Per-pTγ projections + 3×3 tables (shape)
              auto Make3x3Table_DphiShape =
                [&](TH2* h2, const string& outPng, const string& titlePrefix)
              {
                if (!h2) return;

                const int nPt = h2->GetXaxis()->GetNbins();
                const int perPage = kNPtBins;
                int page = 0;

                for (int start = 1; start <= nPt; start += perPage)
                {
                  ++page;

                  TCanvas c(
                    TString::Format("c_tbl_absDphi_%s_%s_p%d", ds.label.c_str(), rKey.c_str(), page).Data(),
                    "c_tbl_absDphi", 1500, 1200
                  );
                  c.Divide(3,3, 0.001, 0.001);

                  vector<TH1*> keep;
                  keep.reserve(perPage);

                  for (int k = 0; k < perPage; ++k)
                  {
                    const int ib = start + k;
                    c.cd(k+1);

                    gPad->SetLeftMargin(0.14);
                    gPad->SetRightMargin(0.05);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetTopMargin(0.10);

                    if (ib > nPt)
                    {
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.06);
                      t.DrawLatex(0.20, 0.55, "EMPTY");
                      continue;
                    }

                    TH1D* p = ProjY_AtXbin(h2, ib,
                      TString::Format("p_absDphi_tbl_%s_%s_%d", ds.label.c_str(), rKey.c_str(), ib).Data()
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

                    NormalizeVisible(p);
                    p->SetTitle("");
                    p->SetLineWidth(2);
                    p->SetMarkerStyle(20);
                    p->SetMarkerSize(1.0);
                    p->GetXaxis()->SetTitle("|#Delta#phi| [rad]");
                    p->GetYaxis()->SetTitle("A.U.");
                    p->Draw("E1");

                    const string ptLab = AxisBinLabel(h2->GetXaxis(), ib, "GeV", 0);

                    TLatex tt;
                    tt.SetNDC(true);
                    tt.SetTextFont(42);
                    tt.SetTextAlign(22);
                    tt.SetTextSize(0.060);
                    tt.DrawLatex(0.50, 0.95,
                      TString::Format("%s, p_{T}^{#gamma} = %s", titlePrefix.c_str(), ptLab.c_str()).Data()
                    );

                    keep.push_back(p);
                  }

                  string nameOut;
                  if (nPt <= perPage)
                  {
                    nameOut = outPng;
                  }
                  else
                  {
                    const size_t pos = outPng.find(".png");
                    const string stem = (pos == string::npos) ? outPng : outPng.substr(0, pos);
                    nameOut = TString::Format("%s_page%d.png", stem.c_str(), page).Data();
                  }

                  SaveCanvas(c, JoinPath(dphiProjDir, nameOut));
                  for (auto* h : keep) delete h;
                }
              };

                // 3×3 tables + per-pT bin plots (counts + shape)
                if (H.hRecoDphi)
                {
                  Make3x3Table_DphiShape(H.hRecoDphi,
                    "table3x3_absDphi_RECO_shape.png",
                    "Reco-level |#Delta#phi(#gamma,jet)|"
                  );

                  const int nPtUse = (int)PtBins().size();
                  TAxis* ax = H.hRecoDphi->GetXaxis();

                  for (int i = 0; i < nPtUse; ++i)
                  {
                    const PtBin& b = PtBins()[i];

                    // Map canonical bin [b.lo,b.hi) onto the unfolding axis via bin-center lookup,
                    // then verify edges match (prevents silent mislabeling).
                    const double mid  = 0.5 * (b.lo + b.hi);
                    const int    xbin = ax->FindBin(mid);

                    if (xbin < 1 || xbin > ax->GetNbins())
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Δphi RECO: cannot map canonical pT bin "
                           << b.lo << "-" << b.hi
                           << " to unfolding axis (xbin=" << xbin << ") → skip"
                           << ANSI_RESET << "\n";
                      continue;
                    }

                    const double axlo = ax->GetBinLowEdge(xbin);
                    const double axhi = ax->GetBinUpEdge(xbin);

                    if (std::fabs(axlo - (double)b.lo) > 1e-6 || std::fabs(axhi - (double)b.hi) > 1e-6)
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Δphi RECO: canonical pT bin "
                           << b.lo << "-" << b.hi
                           << " maps to axis bin " << std::fixed << std::setprecision(0) << axlo << "-" << axhi
                           << " (xbin=" << xbin << ") → skip"
                           << ANSI_RESET << "\n";
                      continue;
                    }

                    TH1D* p = ProjY_AtXbin(H.hRecoDphi, xbin,
                      TString::Format("p_absDphi_reco_%s_%s_%s", ds.label.c_str(), rKey.c_str(), b.folder.c_str()).Data()
                    );
                    if (!p) continue;

                    const string perDir = JoinPath(dphiRecoDir, b.folder);
                    EnsureDir(perDir);

                    vector<string> lines = {
                      "RECO inclusive |#Delta#phi(#gamma,jet)|",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data()
                    };

                    DrawAndSaveTH1_Common(ds, p,
                      JoinPath(perDir, "absDphi_counts.png"),
                      "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1"
                    );

                    TH1* ps = CloneTH1(p,
                      TString::Format("p_absDphi_reco_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                    );
                    NormalizeVisible(ps);

                    DrawAndSaveTH1_Common(ds, ps,
                      JoinPath(perDir, "absDphi_shape.png"),
                      "|#Delta#phi| [rad]", "A.U.", lines, false, false, 0.0, "E1"
                    );

                    delete ps;
                    delete p;
                  }
                }
                if (H.hTruthDphi)
                {
                  Make3x3Table_DphiShape(H.hTruthDphi,
                    "table3x3_absDphi_TRUTH_shape.png",
                    "Truth-level |#Delta#phi(#gamma,jet)|"
                  );

                  const int nPtUse = (int)PtBins().size();
                  TAxis* ax = H.hTruthDphi->GetXaxis();

                  for (int i = 0; i < nPtUse; ++i)
                  {
                    const PtBin& b = PtBins()[i];

                    // Map canonical bin [b.lo,b.hi) onto the unfolding axis via bin-center lookup,
                    // then verify edges match (prevents silent mislabeling).
                    const double mid  = 0.5 * (b.lo + b.hi);
                    const int    xbin = ax->FindBin(mid);

                    if (xbin < 1 || xbin > ax->GetNbins())
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Δphi TRUTH: cannot map canonical pT bin "
                           << b.lo << "-" << b.hi
                           << " to unfolding axis (xbin=" << xbin << ") → skip"
                           << ANSI_RESET << "\n";
                      continue;
                    }

                    const double axlo = ax->GetBinLowEdge(xbin);
                    const double axhi = ax->GetBinUpEdge(xbin);

                    if (std::fabs(axlo - (double)b.lo) > 1e-6 || std::fabs(axhi - (double)b.hi) > 1e-6)
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Δphi TRUTH: canonical pT bin "
                           << b.lo << "-" << b.hi
                           << " maps to axis bin " << std::fixed << std::setprecision(0) << axlo << "-" << axhi
                           << " (xbin=" << xbin << ") → skip"
                           << ANSI_RESET << "\n";
                      continue;
                    }

                    TH1D* p = ProjY_AtXbin(H.hTruthDphi, xbin,
                      TString::Format("p_absDphi_truth_%s_%s_%s", ds.label.c_str(), rKey.c_str(), b.folder.c_str()).Data()
                    );
                    if (!p) continue;

                    const string perDir = JoinPath(dphiTruthDir, b.folder);
                    EnsureDir(perDir);

                    vector<string> lines = {
                      "TRUTH inclusive |#Delta#phi(#gamma,jet)|",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data()
                    };

                    DrawAndSaveTH1_Common(ds, p,
                      JoinPath(perDir, "absDphi_counts.png"),
                      "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1"
                    );

                    TH1* ps = CloneTH1(p,
                      TString::Format("p_absDphi_truth_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data()
                    );
                    NormalizeVisible(ps);

                    DrawAndSaveTH1_Common(ds, ps,
                      JoinPath(perDir, "absDphi_shape.png"),
                      "|#Delta#phi| [rad]", "A.U.", lines, false, false, 0.0, "E1"
                    );

                    delete ps;
                    delete p;
                  }
                }
              };

            // --------------------------------------------------------------------------
            // Part C (SIM-only): response de-flattening into block matrices
            //   - This is your "2D within 2D" visualization (requires hResp + hTruth + hReco)
            // --------------------------------------------------------------------------
            auto RunResponseBlockMatrixQA =
              [&](const string& rKey, double R, const string& rOut, const UnfoldHists& H)
            {
              if (!(H.hResp && H.hTruth && H.hReco)) return;

              const string blockDir = JoinPath(rOut, "ResponseBlockMatrix");
              EnsureDir(blockDir);

              const int nPtTruth = H.hTruth->GetXaxis()->GetNbins();
              const int nPtReco  = H.hReco ->GetXaxis()->GetNbins();
              const int nXJTruth = H.hTruth->GetYaxis()->GetNbins();
              const int nXJReco  = H.hReco ->GetYaxis()->GetNbins();

              const int nGlobTruth = (nPtTruth + 2) * (nXJTruth + 2);
              const int nGlobReco  = (nPtReco  + 2) * (nXJReco  + 2);

              if (H.hResp->GetNbinsX() != nGlobTruth || H.hResp->GetNbinsY() != nGlobReco)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] ResponseBlockMatrix skipped: hResp global-bin dimensions do not match hTruth/hReco.\n"
                     << "  dataset=" << ds.label << "  rKey=" << rKey << "\n"
                     << "  hResp nx,ny = " << H.hResp->GetNbinsX() << "," << H.hResp->GetNbinsY() << "\n"
                     << "  expected   = " << nGlobTruth << "," << nGlobReco
                     << ANSI_RESET << "\n";
                return;
              }

              // (A) mapping file
              vector<string> mapLines;
              mapLines.push_back("globalBinMapping_truthReco.txt");
              mapLines.push_back(string("Dataset: ") + ds.label);
              mapLines.push_back(string("rKey: ") + rKey + TString::Format(" (R=%.1f)", R).Data());
              mapLines.push_back("");
              mapLines.push_back("NOTE: global bin indices are the TH2 internal indices returned by TH2::GetBin/FindBin.");
              mapLines.push_back("      We list ALL bins including ROOT underflow/overflow (bin=0 and bin=nbins+1).");
              mapLines.push_back("");

              mapLines.push_back("[TRUTH] gTruth = hTruth->GetBin(ptBin, xJBin)");
              mapLines.push_back("Columns: gTruth  ptBin  xJBin  pT_range  xJ_range");
              for (int ipt = 0; ipt <= nPtTruth + 1; ++ipt)
              {
                for (int ixj = 0; ixj <= nXJTruth + 1; ++ixj)
                {
                  const int g = H.hTruth->GetBin(ipt, ixj);

                  string ptStr;
                  if (ipt == 0) ptStr = "UNDERFLOW";
                  else if (ipt == nPtTruth + 1) ptStr = "OVERFLOW";
                  else ptStr = AxisBinLabel(H.hTruth->GetXaxis(), ipt, "GeV", 0);

                  string xjStr;
                  if (ixj == 0) xjStr = "UNDERFLOW";
                  else if (ixj == nXJTruth + 1) xjStr = "OVERFLOW";
                  else xjStr = AxisBinLabel(H.hTruth->GetYaxis(), ixj, "", 2);

                  mapLines.push_back(
                    TString::Format("%4d  %2d  %2d  %s  %s",
                                    g, ipt, ixj, ptStr.c_str(), xjStr.c_str()).Data()
                  );
                }
              }

              mapLines.push_back("");
              mapLines.push_back("[RECO] gReco = hReco->GetBin(ptBin, xJBin)");
              mapLines.push_back("Columns: gReco  ptBin  xJBin  pT_range  xJ_range");
              for (int ipt = 0; ipt <= nPtReco + 1; ++ipt)
              {
                for (int ixj = 0; ixj <= nXJReco + 1; ++ixj)
                {
                  const int g = H.hReco->GetBin(ipt, ixj);

                  string ptStr;
                  if (ipt == 0) ptStr = "UNDERFLOW";
                  else if (ipt == nPtReco + 1) ptStr = "OVERFLOW";
                  else ptStr = AxisBinLabel(H.hReco->GetXaxis(), ipt, "GeV", 0);

                  string xjStr;
                  if (ixj == 0) xjStr = "UNDERFLOW";
                  else if (ixj == nXJReco + 1) xjStr = "OVERFLOW";
                  else xjStr = AxisBinLabel(H.hReco->GetYaxis(), ixj, "", 2);

                  mapLines.push_back(
                    TString::Format("%4d  %2d  %2d  %s  %s",
                                    g, ipt, ixj, ptStr.c_str(), xjStr.c_str()).Data()
                  );
                }
              }

              WriteTextFile(JoinPath(blockDir, "globalBinMapping_truthReco.txt"), mapLines);

              // (B) consistent z-range
              double zMax = 0.0;
              double zMinPos = std::numeric_limits<double>::max();

              for (int iptT = 1; iptT <= nPtTruth; ++iptT)
              {
                for (int iptR = 1; iptR <= nPtReco; ++iptR)
                {
                  for (int ixjt = 1; ixjt <= nXJTruth; ++ixjt)
                  {
                    const int gT = H.hTruth->GetBin(iptT, ixjt);
                    const int bx = H.hResp->GetXaxis()->FindBin((double)gT);

                    for (int ixjr = 1; ixjr <= nXJReco; ++ixjr)
                    {
                      const int gR = H.hReco->GetBin(iptR, ixjr);
                      const int by = H.hResp->GetYaxis()->FindBin((double)gR);

                      const double v = H.hResp->GetBinContent(bx, by);
                      if (v > zMax) zMax = v;
                      if (v > 0.0 && v < zMinPos) zMinPos = v;
                    }
                  }
                }
              }
              if (!std::isfinite(zMinPos) || zMinPos == std::numeric_limits<double>::max())
              {
                zMinPos = 1e-6;
              }

              // (C) build one subresponse
              auto MakeSubResponse =
                [&](int iptTruth, int iptReco, const string& name)->TH2F*
              {
                TH2F* hsub = new TH2F(
                  name.c_str(), "",
                  nXJTruth, 0.5, nXJTruth + 0.5,
                  nXJReco,  0.5, nXJReco  + 0.5
                );
                hsub->SetDirectory(nullptr);
                hsub->SetStats(0);
                hsub->GetXaxis()->SetTitle("truth x_{J} bin index");
                hsub->GetYaxis()->SetTitle("reco x_{J} bin index");

                for (int ixjt = 1; ixjt <= nXJTruth; ++ixjt)
                {
                  const int gT = H.hTruth->GetBin(iptTruth, ixjt);
                  const int bx = H.hResp->GetXaxis()->FindBin((double)gT);

                  for (int ixjr = 1; ixjr <= nXJReco; ++ixjr)
                  {
                    const int gR = H.hReco->GetBin(iptReco, ixjr);
                    const int by = H.hResp->GetYaxis()->FindBin((double)gR);

                    hsub->SetBinContent(ixjt, ixjr, H.hResp->GetBinContent(bx, by));
                    hsub->SetBinError  (ixjt, ixjr, H.hResp->GetBinError  (bx, by));
                  }
                }
                return hsub;
              };

              // (D) draw block matrix
              auto DrawBlockMatrix =
                [&](bool logz, const string& outStem)
              {
                const int cw = std::max(1400, 220 * nPtTruth);
                const int ch = std::max(1200, 220 * nPtReco);

                TCanvas c(
                  TString::Format("c_respBlocks_%s_%s_%s",
                                  ds.label.c_str(), rKey.c_str(), logz ? "logz" : "lin").Data(),
                  "c_respBlocks",
                  cw, ch
                );
                c.Divide(nPtTruth, nPtReco, 0.001, 0.001);

                vector<TH2*> keep;
                keep.reserve(nPtTruth * nPtReco);

                for (int iptR = 1; iptR <= nPtReco; ++iptR)
                {
                  for (int iptT = 1; iptT <= nPtTruth; ++iptT)
                  {
                    const int padIdx = (iptR - 1) * nPtTruth + iptT;
                    c.cd(padIdx);

                    gPad->SetLeftMargin(0.12);
                    gPad->SetRightMargin(0.03);
                    gPad->SetBottomMargin(0.12);
                    gPad->SetTopMargin(0.10);
                    gPad->SetTicks(1,1);
                    gPad->SetLogz(logz);

                    TH2F* hsub = MakeSubResponse(
                      iptT, iptR,
                      TString::Format("h_respBlock_%s_%s_t%d_r%d_%s",
                                      ds.label.c_str(), rKey.c_str(), iptT, iptR, logz ? "logz" : "lin").Data()
                    );
                    if (!hsub) continue;

                    if (zMax > 0.0) hsub->SetMaximum(zMax);
                    if (logz) hsub->SetMinimum(std::max(0.5 * zMinPos, 1e-6));
                    else      hsub->SetMinimum(0.0);

                    const bool showX = (iptR == nPtReco);
                    const bool showY = (iptT == 1);

                    hsub->GetXaxis()->SetLabelSize(showX ? 0.08 : 0.0);
                    hsub->GetYaxis()->SetLabelSize(showY ? 0.08 : 0.0);
                    hsub->GetXaxis()->SetTitleSize(showX ? 0.09 : 0.0);
                    hsub->GetYaxis()->SetTitleSize(showY ? 0.09 : 0.0);
                    hsub->GetXaxis()->SetTitleOffset(0.90);
                    hsub->GetYaxis()->SetTitleOffset(0.90);

                    hsub->SetTitle("");
                    hsub->Draw("COL");

                    if (iptR == 1)
                    {
                      const string ptLab = AxisBinLabel(H.hTruth->GetXaxis(), iptT, "GeV", 0);
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.09);
                      t.DrawLatex(0.05, 0.92,
                        TString::Format("T%d %s", iptT, ptLab.c_str()).Data()
                      );
                    }

                    if (iptT == 1)
                    {
                      const string ptLab = AxisBinLabel(H.hReco->GetXaxis(), iptR, "GeV", 0);
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.09);
                      t.DrawLatex(0.05, 0.82,
                        TString::Format("R%d %s", iptR, ptLab.c_str()).Data()
                      );
                    }

                    keep.push_back(hsub);
                  }
                }

                SaveCanvas(c, JoinPath(blockDir, outStem + (logz ? "_logz.png" : "_lin.png")));
                SaveCanvas(c, JoinPath(blockDir, outStem + (logz ? "_logz.pdf" : "_lin.pdf")));

                for (auto* h : keep) delete h;
              };

              DrawBlockMatrix(false, "responseBlockMatrix_xJreco_vs_xJtruth_byPtBins");
              DrawBlockMatrix(true,  "responseBlockMatrix_xJreco_vs_xJtruth_byPtBins");
            };

            // --------------------------------------------------------------------------
            // Part D (SIM only in practice, but safe on data): efficiency/purity summary,
            // 2D maps, xJ distributions, overlays, and summary text.
            //   - This is your big terminal table + all outputs under efficiencyPurity/
            // --------------------------------------------------------------------------
            auto RunEfficiencyPurityQA =
              [&](const string& rKey, double R, const string& rOut, const UnfoldHists& H)
            {
              if (!(H.hReco && H.hTruth)) return;

                const int nPtReco  = H.hReco->GetXaxis()->GetNbins();
                const int nPtTruth = H.hTruth->GetXaxis()->GetNbins();
                const int nPtCanon = (int)PtBins().size();

                // Map canonical pT bin [lo,hi) to the unfolding axis bin by bin-center lookup,
                // then verify the axis bin edges match exactly.
                auto MapCanonicalToAxisBin =
                  [&](TAxis* ax, const PtBin& b, const string& who)->int
                {
                  if (!ax) return -1;

                  const double mid  = 0.5 * (b.lo + b.hi);
                  const int    xbin = ax->FindBin(mid);

                  if (xbin < 1 || xbin > ax->GetNbins()) return -1;

                  const double axlo = ax->GetBinLowEdge(xbin);
                  const double axhi = ax->GetBinUpEdge(xbin);

                  if (std::fabs(axlo - (double)b.lo) > 1e-6 || std::fabs(axhi - (double)b.hi) > 1e-6)
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] Unfolding pT mapping mismatch (" << who << "): requested "
                         << b.lo << "-" << b.hi
                         << " but axis bin " << xbin << " is "
                         << std::fixed << std::setprecision(0) << axlo << "-" << axhi
                         << " → skipping this pT bin"
                         << ANSI_RESET << "\n";
                    return -1;
                  }
                  return xbin;
                };

                const string dirEP = JoinPath(rOut, "efficiencyPurity");
                const string dirEP_Maps = JoinPath(dirEP, "maps2D");
                const string dirEP_XJ   = JoinPath(dirEP, "xJ_distributions");
                EnsureDir(dirEP);
                EnsureDir(dirEP_Maps);
                EnsureDir(dirEP_XJ);

                cout << ANSI_BOLD_CYN << "\n[Unfolding table] " << ds.label << "  rKey=" << rKey
                     << " (R=" << std::fixed << std::setprecision(1) << R << ")\n" << ANSI_RESET;

                const int wPt = 16;
                const int wN  = 14;
                const int wF  = 12;

                cout << std::left << std::setw(wPt) << "pTgamma bin"
                     << std::right
                     << std::setw(wN) << "N_truth"
                     << std::setw(wN) << "N_miss"
                     << std::setw(wF) << "eff"
                     << std::setw(wN) << "N_reco"
                     << std::setw(wN) << "N_fake"
                     << std::setw(wF) << "pur"
                     << "\n";
                cout << string(wPt + 4*wN + 2*wF, '-') << "\n";

                vector<string> lines;
                lines.push_back(string("Unfolding efficiency/purity summary (") + ds.label + ")");
                lines.push_back(string("rKey: ") + rKey + TString::Format("  R=%.1f", R).Data());
                lines.push_back("");
                lines.push_back("Definitions:");
                lines.push_back("  efficiency(truth->reco) = (N_truth - N_miss) / N_truth");
                lines.push_back("  purity(reco sample)     = (N_reco  - N_fake) / N_reco");
                lines.push_back("");

                vector<double> xCenters;
                vector<double> effVsPt;
                vector<double> purVsPt;
                xCenters.reserve(nPtCanon);
                effVsPt.reserve(nPtCanon);
                purVsPt.reserve(nPtCanon);

              // 2D efficiency map (if misses exist)
              if (H.hMisses)
              {
                TH2* hEff2D = CloneTH2(H.hTruth, "h2_efficiency_pTgamma_xJ");
                if (hEff2D)
                {
                  hEff2D->Reset("ICES");
                  const int nx = H.hTruth->GetXaxis()->GetNbins();
                  const int ny = H.hTruth->GetYaxis()->GetNbins();

                  for (int ix = 1; ix <= nx; ++ix)
                  {
                    for (int iy = 1; iy <= ny; ++iy)
                    {
                      const double t = H.hTruth ->GetBinContent(ix, iy);
                      const double m = H.hMisses->GetBinContent(ix, iy);
                      const double eff = (t > 0.0) ? ((t - m) / t) : 0.0;
                      hEff2D->SetBinContent(ix, iy, eff);
                    }
                  }

                  hEff2D->SetMinimum(0.0);
                  hEff2D->SetMaximum(1.0);

                  DrawAndSaveTH2_Common(ds, hEff2D,
                    JoinPath(dirEP_Maps, "efficiency2D_pTgamma_vs_xJ.png"),
                    "p_{T}^{#gamma,truth} [GeV]", "x_{J#gamma}^{truth}", "Efficiency",
                    { "Efficiency map from unfolding inputs",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      "eff = (truth - misses) / truth" },
                    false);

                  delete hEff2D;
                }
              }

              // 2D purity map (if fakes exist)
              if (H.hFakes)
              {
                TH2* hPur2D = CloneTH2(H.hReco, "h2_purity_pTgamma_xJ");
                if (hPur2D)
                {
                  hPur2D->Reset("ICES");
                  const int nx = H.hReco->GetXaxis()->GetNbins();
                  const int ny = H.hReco->GetYaxis()->GetNbins();

                  for (int ix = 1; ix <= nx; ++ix)
                  {
                    for (int iy = 1; iy <= ny; ++iy)
                    {
                      const double r  = H.hReco ->GetBinContent(ix, iy);
                      const double fk = H.hFakes->GetBinContent(ix, iy);
                      const double pur = (r > 0.0) ? ((r - fk) / r) : 0.0;
                      hPur2D->SetBinContent(ix, iy, pur);
                    }
                  }

                  hPur2D->SetMinimum(0.0);
                  hPur2D->SetMaximum(1.0);

                  DrawAndSaveTH2_Common(ds, hPur2D,
                    JoinPath(dirEP_Maps, "purity2D_pTgamma_vs_xJ.png"),
                    "p_{T}^{#gamma,reco} [GeV]", "x_{J#gamma}^{reco}", "Purity",
                    { "Purity map from unfolding inputs",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      "pur = (reco - fakes) / reco" },
                    false);

                  delete hPur2D;
                }
              }

              // 3x3 tables metric vs xJ
              auto Make3x3Table_MetricVsXJ =
                [&](TH2* hDen, TH2* hNum,
                    const string& outPng,
                    const string& metricTitle,
                    const string& metricExpr,
                    bool isTruthMetric)
              {
                if (!hDen || !hNum) return;

                const int nx = hDen->GetXaxis()->GetNbins();
                const int perPage = 9;

                int page = 0;
                for (int start = 1; start <= nx; start += perPage)
                {
                  ++page;

                  TCanvas c(
                    TString::Format("c_tbl_%s_%s_p%d", metricTitle.c_str(), rKey.c_str(), page).Data(),
                    "c_tbl_metric", 1500, 1200
                  );
                  c.Divide(3,3, 0.001, 0.001);

                  std::vector<TH1*> keep;
                  keep.reserve(perPage);

                  for (int k = 0; k < perPage; ++k)
                  {
                    const int ix = start + k;
                    c.cd(k+1);

                    gPad->SetLeftMargin(0.14);
                    gPad->SetRightMargin(0.05);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetTopMargin(0.10);
                    gPad->SetLogy(false);

                    if (ix > nx)
                    {
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.06);
                      t.DrawLatex(0.20, 0.55, "EMPTY");
                      continue;
                    }

                    TH1D* hDenY = hDen->ProjectionY(
                      TString::Format("den_%s_%d", metricTitle.c_str(), ix).Data(), ix, ix
                    );
                    TH1D* hNumY = hNum->ProjectionY(
                      TString::Format("num_%s_%d", metricTitle.c_str(), ix).Data(), ix, ix
                    );

                    if (!hDenY || !hNumY)
                    {
                      if (hDenY) delete hDenY;
                      if (hNumY) delete hNumY;
                      TLatex t;
                      t.SetNDC(true);
                      t.SetTextFont(42);
                      t.SetTextSize(0.06);
                      t.DrawLatex(0.15, 0.55, "MISSING");
                      continue;
                    }

                    hDenY->SetDirectory(nullptr);
                    hNumY->SetDirectory(nullptr);

                    TH1D* hMet = (TH1D*)hDenY->Clone(TString::Format("met_%s_%d", metricTitle.c_str(), ix).Data());
                    hMet->SetDirectory(nullptr);
                    hMet->Reset("ICES");

                    const int ny = hMet->GetNbinsX();
                    for (int iy = 1; iy <= ny; ++iy)
                    {
                      const double den = hDenY->GetBinContent(iy);
                      const double num = hNumY->GetBinContent(iy);
                      const double val = (den > 0.0) ? ((den - num) / den) : 0.0;
                      const double err = (den > 0.0) ? std::sqrt(std::max(0.0, val*(1.0 - val)/den)) : 0.0;

                      hMet->SetBinContent(iy, val);
                      hMet->SetBinError(iy, err);
                    }

                    hMet->SetMinimum(0.0);
                    hMet->SetMaximum(1.05);
                    hMet->SetLineWidth(2);
                    hMet->SetMarkerStyle(20);
                    hMet->SetTitle("");
                    hMet->GetXaxis()->SetTitle(isTruthMetric ? "x_{J#gamma}^{truth}" : "x_{J#gamma}^{reco}");
                    hMet->GetYaxis()->SetTitle(metricTitle.c_str());
                    hMet->Draw("E1");

                    const string ptLab = AxisBinLabel(hDen->GetXaxis(), ix, "GeV", 0);

                    vector<string> box;
                    box.push_back(metricExpr);
                    box.push_back(TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data());
                    box.push_back(TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());
                    DrawLatexLines(0.16, 0.90, box, 0.040, 0.050);

                    keep.push_back(hMet);

                    delete hDenY;
                    delete hNumY;
                  }

                  string nameOut;
                  if (nx <= perPage)
                  {
                    nameOut = outPng;
                  }
                  else
                  {
                    const size_t pos = outPng.find(".png");
                    const string stem = (pos == string::npos) ? outPng : outPng.substr(0, pos);
                    nameOut = TString::Format("%s_page%d.png", stem.c_str(), page).Data();
                  }

                  SaveCanvas(c, JoinPath(dirEP_XJ, nameOut));
                  for (auto* h : keep) delete h;
                }
              };

              if (H.hMisses)
              {
                Make3x3Table_MetricVsXJ(
                  H.hTruth, H.hMisses,
                  "table3x3_efficiency_vs_xJ.png",
                  "Efficiency",
                  "eff(x_{J}) = (truth - misses)/truth",
                  true
                );
              }
              if (H.hFakes)
              {
                Make3x3Table_MetricVsXJ(
                  H.hReco, H.hFakes,
                  "table3x3_purity_vs_xJ.png",
                  "Purity",
                  "pur(x_{J}) = (reco - fakes)/reco",
                  false
                );
              }

                // Per canonical pT bin: terminal line + text summary + overlay truth vs reco xJ shapes
                for (int i = 0; i < nPtCanon; ++i)
                {
                  const PtBin& b = PtBins()[i];

                  const int xTruth = MapCanonicalToAxisBin(H.hTruth->GetXaxis(), b, "truth");
                  const int xReco  = MapCanonicalToAxisBin(H.hReco ->GetXaxis(), b, "reco");
                  if (xTruth < 1 || xReco < 1) continue;

                  const double xC = 0.5 * (b.lo + b.hi);
                  xCenters.push_back(xC);

                  const string ptLabTruth = AxisBinLabel(H.hTruth->GetXaxis(), xTruth, "GeV", 0);
                  const string ptLabReco  = AxisBinLabel(H.hReco ->GetXaxis(), xReco,  "GeV", 0);

                  const double Ntruth = H.hTruth->Integral(xTruth, xTruth, 0, H.hTruth->GetYaxis()->GetNbins() + 1);
                  const double Nreco  = H.hReco ->Integral(xReco,  xReco,  0, H.hReco ->GetYaxis()->GetNbins() + 1);

                  const double Nmiss  = (H.hMisses ? H.hMisses->Integral(xTruth, xTruth, 0, H.hMisses->GetYaxis()->GetNbins() + 1) : 0.0);
                  const double Nfake  = (H.hFakes  ? H.hFakes ->Integral(xReco,  xReco,  0, H.hFakes ->GetYaxis()->GetNbins() + 1) : 0.0);

                  const double eff = SafeDivide(Ntruth - Nmiss, Ntruth, 0.0);
                  const double pur = SafeDivide(Nreco  - Nfake, Nreco,  0.0);

                  effVsPt.push_back(eff);
                  purVsPt.push_back(pur);

                  cout << std::left << std::setw(wPt) << ptLabTruth
                       << std::right
                       << std::setw(wN) << std::fixed << std::setprecision(0) << Ntruth
                       << std::setw(wN) << Nmiss
                       << std::setw(wF) << std::fixed << std::setprecision(4) << eff
                       << std::setw(wN) << std::fixed << std::setprecision(0) << Nreco
                       << std::setw(wN) << Nfake
                       << std::setw(wF) << std::fixed << std::setprecision(4) << pur
                       << "\n";

                  lines.push_back(TString::Format(
                    "pTtruth=%s  pTreco=%s  Ntruth=%.0f  Nmiss=%.0f  eff=%.6f  Nreco=%.0f  Nfake=%.0f  pur=%.6f",
                    ptLabTruth.c_str(), ptLabReco.c_str(),
                    Ntruth, Nmiss, eff,
                    Nreco, Nfake, pur
                  ).Data());

                  // Truth-vs-reco xJ shape overlay per canonical pT bin (aligned by edges)
                  {
                    TH1D* pxTruth = H.hTruth->ProjectionY(
                      TString::Format("pxTruth_%s_%d", rKey.c_str(), i).Data(), xTruth, xTruth
                    );
                    TH1D* pxReco  = H.hReco ->ProjectionY(
                      TString::Format("pxReco_%s_%d", rKey.c_str(), i).Data(),  xReco,  xReco
                    );

                    if (pxTruth && pxReco)
                    {
                      pxTruth->SetDirectory(nullptr);
                      pxReco->SetDirectory(nullptr);
                      NormalizeToUnitArea(pxTruth);
                      NormalizeToUnitArea(pxReco);

                      pxTruth->SetLineWidth(2);
                      pxReco->SetLineWidth(2);
                      pxTruth->SetLineColor(1);
                      pxReco->SetLineColor(2);

                      TCanvas c(
                        TString::Format("c_unf_ov_%s_%d", rKey.c_str(), i).Data(),
                        "c_unf_ov", 900, 700
                      );
                      ApplyCanvasMargins1D(c);

                      const double maxv = std::max(pxTruth->GetMaximum(), pxReco->GetMaximum());
                      pxTruth->SetMaximum(maxv * 1.25);
                      pxTruth->SetTitle("");
                      pxTruth->GetXaxis()->SetTitle("x_{J}");
                      pxTruth->GetYaxis()->SetTitle("A.U.");
                      pxTruth->Draw("hist");
                      pxReco->Draw("hist same");

                      TLegend leg(0.60, 0.78, 0.92, 0.90);
                      leg.SetTextFont(42);
                      leg.SetTextSize(0.033);
                      leg.AddEntry(pxTruth, "Truth (shape)", "l");
                      leg.AddEntry(pxReco,  "Reco (shape)",  "l");
                      leg.Draw();

                      DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                      DrawLatexLines(
                        0.14, 0.78,
                        { string("Truth vs reco x_{J} shape"), rKey, TString::Format("p_{T}^{#gamma}: %s", ptLabTruth.c_str()).Data() },
                        0.030, 0.040
                      );

                      SaveCanvas(c, JoinPath(rOut, TString::Format("overlay_truth_vs_reco_xJ_shape_pTbin%d.png", i + 1).Data()));

                      delete pxTruth;
                      delete pxReco;
                    }
                    else
                    {
                      if (pxTruth) delete pxTruth;
                      if (pxReco)  delete pxReco;
                    }
                  }
                }

                // Presentation graphs: efficiency and purity vs pTgamma (integrated over xJ)
                {
                  const int nPtPlot = (int)xCenters.size();
                  if (nPtPlot > 0)
                  {
                    TCanvas c1(TString::Format("c_eff_vs_pt_%s", rKey.c_str()).Data(), "c_eff_vs_pt", 900,700);
                    ApplyCanvasMargins1D(c1);
                    TGraph g(nPtPlot, &xCenters[0], &effVsPt[0]);
                    g.SetLineWidth(2);
                    g.SetMarkerStyle(20);
                    g.Draw("ALP");
                    g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                    g.GetYaxis()->SetTitle("Efficiency (integrated over x_{J})");
                    g.SetMinimum(0.0);
                    g.SetMaximum(1.05);

                    DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                    DrawLatexLines(0.14,0.78, { "Unfolding efficiency vs p_{T}^{#gamma}", rKey }, 0.030, 0.040);

                    SaveCanvas(c1, JoinPath(dirEP, "efficiency_vs_pTgamma_integratedXJ.png"));
                  }
                }
                {
                  const int nPtPlot = (int)xCenters.size();
                  if (nPtPlot > 0)
                  {
                    TCanvas c2(TString::Format("c_pur_vs_pt_%s", rKey.c_str()).Data(), "c_pur_vs_pt", 900,700);
                    ApplyCanvasMargins1D(c2);
                    TGraph g(nPtPlot, &xCenters[0], &purVsPt[0]);
                    g.SetLineWidth(2);
                    g.SetMarkerStyle(20);
                    g.Draw("ALP");
                    g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                    g.GetYaxis()->SetTitle("Purity (integrated over x_{J})");
                    g.SetMinimum(0.0);
                    g.SetMaximum(1.05);

                    DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                    DrawLatexLines(0.14,0.78, { "Unfolding purity vs p_{T}^{#gamma}", rKey }, 0.030, 0.040);

                    SaveCanvas(c2, JoinPath(dirEP, "purity_vs_pTgamma_integratedXJ.png"));
                  }
                }


              WriteTextFile(JoinPath(rOut, "summary_unfolding_eff_pur.txt"), lines);
            };

            // --------------------------------------------------------------------------
            // MAIN LOOP over rKeys (orchestrator)
            // --------------------------------------------------------------------------
            for (const auto& rKey : kRKeys)
            {
              const double R = RFromKey(rKey);

              const string rOut = JoinPath(outDir, rKey);
              EnsureDir(rOut);

              UnfoldHists H = LoadUnfoldHists(rKey);

              if (!HasAnyUnfold(H))
              {
                cout << ANSI_BOLD_YEL << "[WARN] No unfolding histograms found for " << ds.label << " " << rKey << ANSI_RESET << "\n";
                continue;
              }

              // Part A: core unfolding maps + response scatter (SIM + DATA)
              RunCoreUnfolding2D(rKey, R, rOut, H);

              // Part B: Δphi unfolding QA (runs only if histograms exist)
              RunDeltaPhiUnfoldingQA(rKey, R, rOut, H);

              // Part C: response "2D within 2D" de-flattening (requires hResp + hTruth + hReco)
              RunResponseBlockMatrixQA(rKey, R, rOut, H);

              // Part D: efficiency/purity products + overlays + summary file (requires hReco + hTruth)
              RunEfficiencyPurityQA(rKey, R, rOut, H);
           }
        }
        // =============================================================================
        // High-level runner
        // =============================================================================
        // NOTE: Run-mode validation now lives in AnalyzeRecoilJets.h (ValidateRunConfig()).
        // This wrapper is kept only to avoid stale references inside the analysis namespace.
        static bool ExactlyOneModeSet()
        {
            return ValidateRunConfig(nullptr);
        }


  } // namespace analysis


  // =============================================================================
  // Driver
  // =============================================================================
  namespace driver
  {
    // When running multiple SIM selections sequentially in SIM+DATA mode,
    // we need a unique PP output base per SIM selection to avoid overwriting.
    // If empty, PP outputs go to kOutPPBase (legacy behavior).
    static string gPPOutBaseSubdir = "";

    inline void SetPPOutBaseSubdir(const string& subdir)
    {
        gPPOutBaseSubdir = subdir;
    }

    inline string PPOutBaseForThisRun()
    {
        if (gPPOutBaseSubdir.empty()) return kOutPPBase;
        return JoinPath(kOutPPBase, gPPOutBaseSubdir);
    }

    inline bool OpenDataset(Dataset& ds)
    {
        // ------------------------------------------------------------------
        // Guard: never allow relative-path output (prevents creating ./RecoilJetQA
        // under whatever directory you launched ROOT from, e.g. macros/).
        // ------------------------------------------------------------------
        if (ds.outBase.empty())
        {
          cout << ANSI_BOLD_RED
               << "[FATAL] Dataset '" << ds.label << "' has EMPTY outBase; refusing relative output.\n"
               << "        inFilePath = " << ds.inFilePath
               << ANSI_RESET << "\n";
          return false;
        }
        if (ds.outBase[0] != '/')
        {
          cout << ANSI_BOLD_RED
               << "[FATAL] Dataset '" << ds.label << "' has NON-ABSOLUTE outBase: " << ds.outBase << "\n"
               << "        Refusing to write outputs relative to the current working directory.\n"
               << "        inFilePath = " << ds.inFilePath
               << ANSI_RESET << "\n";
          return false;
        }

        EnsureDir(ds.outBase);

        ds.file = TFile::Open(ds.inFilePath.c_str(), "READ");
        if (!ds.file || ds.file->IsZombie())
        {
          cout << ANSI_BOLD_RED << "[FATAL] Cannot open input file: " << ds.inFilePath
               << ANSI_RESET << "\n";
          return false;
        }

        ds.topDir = ds.file->GetDirectory(ds.topDirName.c_str());
        if (!ds.topDir)
        {
          cout << ANSI_BOLD_RED << "[FATAL] Missing topDir '" << ds.topDirName
               << "' in file: " << ds.inFilePath << ANSI_RESET << "\n";
          return false;
        }

        const string missPath = JoinPath(ds.outBase, "missing_hists_" + ds.label + ".txt");
        ds.missingOut.open(missPath.c_str());
        ds.missingCount = 0;

        // reset coverage tracking
        ds.requestCounts.clear();
        ds.missingCounts.clear();
        ds.missingReason.clear();

        return true;
    }

    inline void CloseDataset(Dataset& ds)
    {
        if (ds.missingOut.is_open()) ds.missingOut.close();
        if (ds.file) { ds.file->Close(); ds.file = nullptr; ds.topDir = nullptr; }
    }

    inline vector<Dataset> BuildDatasets(RunMode mode)
    {
        vector<Dataset> datasets;

        const bool doSim = (mode == RunMode::kSimOnly ||
                            mode == RunMode::kSimAndDataPP);

        const bool doPP  = (mode == RunMode::kPPDataOnly ||
                            mode == RunMode::kSimAndDataPP);

        const SimSample ss = CurrentSimSample();

        if (doSim)
        {
          Dataset ds;
          ds.label      = "SIM";
          ds.isSim      = true;
          ds.trigger    = "";
          ds.topDirName = kDirSIM;
          ds.inFilePath = SimInputPathForSample(ss);
          ds.outBase    = SimOutBaseForSample(ss);
          datasets.push_back(std::move(ds));
        }

        if (doPP)
        {
          Dataset ds;
          ds.label      = "DATA_PP";
          ds.isSim      = false;
          ds.trigger    = kTriggerPP;
          ds.topDirName = kTriggerPP;
          ds.inFilePath = kInPP;

            // NOTE:
            // - Legacy single-run behavior: PP output goes to kOutPPBase
            // - Multi-run SIM+DATA behavior: PP output goes to kOutPPBase/with_<simSampleLabel>
            //   (so you can loop over multiple SIM selections without clobbering PP outputs)
            //
            // DATA outputs are organized per-trigger:
            //   <PP base>/<trigger>/{baselineData,insituCalib,unfolding,...}
            if (mode == RunMode::kPPDataOnly)
            {
              ds.outBase = JoinPath(kOutPPBase, ds.trigger);
            }
            else
            {
              ds.outBase = JoinPath(PPOutBaseForThisRun(), ds.trigger);
            }

          datasets.push_back(std::move(ds));
        }

        return datasets;
    }

    inline bool MaybeBuildMergedSIM(RunMode mode)
    {
      // Only relevant when a SIM-including mode is running AND a merged SIM sample was selected.
      if (mode == RunMode::kPPDataOnly) return true;

      const SimSample ss = CurrentSimSample();
      if (!IsMergedSimSample(ss)) return true;

      bool ok = true;

      if (ss == SimSample::kPhotonJet5And10Merged)
      {
          ok = BuildMergedSIMFile_PhotonSlices(
            {kInSIM5, DefaultSim10and20Config().photon10},
            {kSigmaPhoton5_pb, kSigmaPhoton10_pb},
            kMergedSIMOut_5and10,
            kDirSIM,
            {"photonJet5", "photonJet10"}
          );
      }
      else if (ss == SimSample::kPhotonJet5And20Merged)
      {
          ok = BuildMergedSIMFile_PhotonSlices(
            {kInSIM5, DefaultSim10and20Config().photon20},
            {kSigmaPhoton5_pb, kSigmaPhoton20_pb},
            kMergedSIMOut_5and20,
            kDirSIM,
            {"photonJet5", "photonJet20"}
          );
      }
      else if (ss == SimSample::kPhotonJet10And20Merged)
      {
          if (!doRemergePhoton10and20sim)
          {
              const string outMerged = MergedSIMOut_10and20_Default();

              cout << ANSI_BOLD_CYN
                   << "\n[MERGE SIM] doRemergePhoton10and20sim=false -> skipping SIM10+20 rebuild step.\n"
                   << "            Using existing merged output: " << outMerged << "\n"
                   << ANSI_RESET;

              if (gSystem->AccessPathName(outMerged.c_str()))
              {
                cout << ANSI_BOLD_RED
                     << "[MERGE SIM][FATAL] Merged SIM10+20 file not found, but doRemergePhoton10and20sim=false:\n"
                     << "  " << outMerged << "\n"
                     << ANSI_RESET;
                ok = false;
              }
          }
          else
          {
              const auto& cfgs = Sim10and20Configs();
              for (const auto& kv : cfgs)
              {
                  const auto& cfg = kv.second;

                  const string outMerged =
                      MergedSIMOut_10and20_ForKey(cfg.key);

                  ok = BuildMergedSIMFile_PhotonSlices(
                    {cfg.photon10, cfg.photon20},
                    {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                    outMerged,
                    kDirSIM,
                    {"photonJet10", "photonJet20"}
                  );

                  if (!ok) break;
              }
          }
      }
      else if (ss == SimSample::kPhotonJet5And10And20Merged)
      {
          ok = BuildMergedSIMFile_PhotonSlices(
            {kInSIM5, DefaultSim10and20Config().photon10, DefaultSim10and20Config().photon20},
            {kSigmaPhoton5_pb, kSigmaPhoton10_pb, kSigmaPhoton20_pb},
            kMergedSIMOut_5and10and20,
            kDirSIM,
            {"photonJet5", "photonJet10", "photonJet20"}
          );
      }

      if (!ok)
      {
        cout << ANSI_BOLD_RED << "[FATAL] Failed to build merged SIM file." << ANSI_RESET << "\n";
        return false;
      }

      return true;
    }
  } // namespace driver


  // =============================================================================
  // ARJ::Run()  (ROOT macro entrypoint called by global AnalyzeRecoilJets())
  // =============================================================================
  int Run()
  {
      SetupGlobalStyle();

      // ---------------------------------------------------------------------------
      // Banner
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN
           << "\n============================================================\n"
           << " AnalyzeRecoilJets  (refactored driver + modular analysis)\n"
           << "============================================================\n"
           << ANSI_RESET;

      // ---------------------------------------------------------------------------
      // Mode validation + printout
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 0] Run-mode validation\n" << ANSI_RESET;

      cout << "  Toggles:\n"
           << "    isPPdataOnly            = " << (isPPdataOnly ? "true" : "false") << "\n"
           << "    isSimAndDataPP          = " << (isSimAndDataPP ? "true" : "false") << "\n"
           << "    isPhotonJet5            = " << (isPhotonJet5 ? "true" : "false") << "\n"
           << "    isPhotonJet10           = " << (isPhotonJet10 ? "true" : "false") << "\n"
           << "    isPhotonJet20           = " << (isPhotonJet20 ? "true" : "false") << "\n"
           << "    bothPhoton5and10sim     = " << (bothPhoton5and10sim ? "true" : "false") << "\n"
           << "    bothPhoton5and20sim     = " << (bothPhoton5and20sim ? "true" : "false") << "\n"
           << "    bothPhoton10and20sim    = " << (bothPhoton10and20sim ? "true" : "false") << "\n"
           << "    doRemergePhoton10and20sim = " << (doRemergePhoton10and20sim ? "true" : "false") << "\n"
           << "    allPhoton5and10and20sim = " << (allPhoton5and10and20sim ? "true" : "false") << "\n";

      // ---------------------------------------------------------------------------
      //  Multi-SIM sequential runner
      //   If multiple SIM sample toggles are true, we automatically run each selected
      //   SIM sample one-after-another (SIM_ONLY or SIM+DATA_PP).
      //   Each run is executed with exactly ONE SIM sample enabled so the rest of the
      //   pipeline remains unchanged.
      // ---------------------------------------------------------------------------
      auto RequestedSamplesFromFlags = [&]() -> vector<SimSample>
      {
        vector<SimSample> v;

        // Order matters (matches the way you'd normally run these by hand):
        //   singles (5,10,20) -> pair merges -> 3-way merge
        if (isPhotonJet5)            v.push_back(SimSample::kPhotonJet5);
        if (isPhotonJet10)           v.push_back(SimSample::kPhotonJet10);
        if (isPhotonJet20)           v.push_back(SimSample::kPhotonJet20);
        if (bothPhoton5and10sim)     v.push_back(SimSample::kPhotonJet5And10Merged);
        if (bothPhoton5and20sim)     v.push_back(SimSample::kPhotonJet5And20Merged);
        if (bothPhoton10and20sim)    v.push_back(SimSample::kPhotonJet10And20Merged);
        if (allPhoton5and10and20sim) v.push_back(SimSample::kPhotonJet5And10And20Merged);

        return v;
      };

      static bool sInMultiRun = false;

      const vector<SimSample> requested = RequestedSamplesFromFlags();

      if (!sInMultiRun)
      {
        // Clear any stale override from a previous call in the same ROOT session.
        driver::SetPPOutBaseSubdir("");

        // Multi-run is only meaningful for SIM-including modes.
        if (!isPPdataOnly && requested.size() > 1)
        {
          cout << ANSI_BOLD_CYN
               << "\n[INFO] Multiple SIM sample toggles are TRUE (N=" << requested.size() << ").\n"
               << "       Will run each selected SIM sample sequentially in mode: "
               << (isSimAndDataPP ? "SIM_AND_DATA_PP" : "SIM_ONLY")
               << "\n"
               << ANSI_RESET;

          // Save original flags so we can restore at the end.
          const bool o_isPhotonJet5            = isPhotonJet5;
          const bool o_isPhotonJet10           = isPhotonJet10;
          const bool o_isPhotonJet20           = isPhotonJet20;
          const bool o_bothPhoton5and10sim     = bothPhoton5and10sim;
          const bool o_bothPhoton5and20sim     = bothPhoton5and20sim;
          const bool o_bothPhoton10and20sim    = bothPhoton10and20sim;
          const bool o_allPhoton5and10and20sim = allPhoton5and10and20sim;

          auto ClearAllSimFlags = [&]()
          {
            isPhotonJet5            = false;
            isPhotonJet10           = false;
            isPhotonJet20           = false;
            bothPhoton5and10sim     = false;
            bothPhoton5and20sim     = false;
            bothPhoton10and20sim    = false;
            allPhoton5and10and20sim = false;
          };

          auto SetFlagsForSample = [&](SimSample s)
          {
            ClearAllSimFlags();
            switch (s)
            {
              case SimSample::kPhotonJet5:                 isPhotonJet5 = true; break;
              case SimSample::kPhotonJet10:                isPhotonJet10 = true; break;
              case SimSample::kPhotonJet20:                isPhotonJet20 = true; break;
              case SimSample::kPhotonJet5And10Merged:      bothPhoton5and10sim = true; break;
              case SimSample::kPhotonJet5And20Merged:      bothPhoton5and20sim = true; break;
              case SimSample::kPhotonJet10And20Merged:     bothPhoton10and20sim = true; break;
              case SimSample::kPhotonJet5And10And20Merged: allPhoton5and10and20sim = true; break;
              default: break;
            }
          };

          int worstRc = 0;

          sInMultiRun = true;
          for (size_t i = 0; i < requested.size(); ++i)
          {
            const SimSample s = requested[i];

            cout << ANSI_BOLD_CYN
                 << "\n============================================================\n"
                 << "[MULTI-RUN " << (i + 1) << "/" << requested.size() << "] SIM sample: " << SimSampleLabel(s) << "\n"
                 << "============================================================\n"
                 << ANSI_RESET;

            SetFlagsForSample(s);

            // If we are doing SIM+DATA, avoid overwriting PP outputs across runs.
            // PP outputs will go under:
            //   kOutPPBase/with_<SimSampleLabel(s)>/
            if (isSimAndDataPP)
            {
              driver::SetPPOutBaseSubdir(string("with_") + SimSampleLabel(s));
            }
            else
            {
              driver::SetPPOutBaseSubdir("");
            }

            const int rc = Run();  // recursion into the single-sample path
            if (rc != 0) worstRc = rc;
          }

          // Restore original flags
          isPhotonJet5            = o_isPhotonJet5;
          isPhotonJet10           = o_isPhotonJet10;
          isPhotonJet20           = o_isPhotonJet20;
          bothPhoton5and10sim     = o_bothPhoton5and10sim;
          bothPhoton5and20sim     = o_bothPhoton5and20sim;
          bothPhoton10and20sim    = o_bothPhoton10and20sim;
          allPhoton5and10and20sim = o_allPhoton5and10and20sim;

          driver::SetPPOutBaseSubdir("");
          sInMultiRun = false;

          return worstRc;
        }
      }

      // ---------------------------------------------------------------------------
      // Normal single-sample path (legacy behavior)
      // ---------------------------------------------------------------------------
      string cfgErr;
      if (!ValidateRunConfig(&cfgErr))
      {
        cout << ANSI_BOLD_RED
             << "[FATAL] Invalid run configuration:\n"
             << "  " << cfgErr << "\n"
             << ANSI_RESET;
        return 1;
      }

      const RunMode mode   = CurrentRunMode();
      const SimSample ss   = CurrentSimSample();

      cout << ANSI_BOLD_YEL << "  -> Selected mode: " << RunModeLabel(mode);
      if (mode != RunMode::kPPDataOnly)
      {
        cout << "  |  SIM sample: " << SimSampleLabel(ss)
             << "  |  SIM outBase: " << SimOutBaseForSample(ss);
      }
      cout << ANSI_RESET << "\n";

      // ---------------------------------------------------------------------------
      // Explicit default SIM10+20 config printout (so it's obvious what drives the merge + analysis)
      // ---------------------------------------------------------------------------
      if (mode != RunMode::kPPDataOnly)
      {
        const Sim10and20Config& cfgDef = DefaultSim10and20Config();

        cout << ANSI_DIM
             << "\n  [SIM10+20 DEFAULT CONFIG]\n"
             << "    DefaultSimSampleKey()           = " << DefaultSimSampleKey() << "\n"
             << "    cfg.key                         = " << cfgDef.key << "\n"
             << "    cfg.photon10                    = " << cfgDef.photon10 << "\n"
             << "    cfg.photon20                    = " << cfgDef.photon20 << "\n"
             << "    cfg.jetMinPt                    = " << cfgDef.jetMinPt << " GeV\n"
             << "    cfg.bbLabel                     = " << cfgDef.bbLabel << "\n"
             << "    kInSIM10 (alias)                = " << kInSIM10 << "\n"
             << "    kInSIM20 (alias)                = " << kInSIM20 << "\n"
             << "    MergedSIMOut_10and20_Default()  = " << MergedSIMOut_10and20_Default() << "\n"
             << ANSI_RESET;

        if (ss == SimSample::kPhotonJet10And20Merged)
        {
          cout << ANSI_DIM
               << "  [SIM10+20 MERGED DATASET PATH]\n"
               << "    SimInputPathForSample(ss)       = " << SimInputPathForSample(ss) << "\n"
               << "    (This is the file opened in STEP 3 and used for the full analysis.)\n"
               << ANSI_RESET;
        }
      }

      // ---------------------------------------------------------------------------
      // Optional SIM slice merge
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 1] Optional SIM slice merge (photonJet5/10/20 combinations)\n" << ANSI_RESET;
      cout << "  -> Checking whether merged SIM is required for this mode...\n";

      if (!driver::MaybeBuildMergedSIM(mode))
      {
        cout << ANSI_BOLD_RED
             << "  [FATAL] driver::MaybeBuildMergedSIM failed. Aborting.\n"
             << ANSI_RESET;
        return 1;
      }
      cout << ANSI_BOLD_GRN << "  [OK] SIM merge step complete (or not needed).\n" << ANSI_RESET;

      // ---------------------------------------------------------------------------
      // Build datasets list
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 2] Build dataset list for this run-mode\n" << ANSI_RESET;

      vector<Dataset> datasets = driver::BuildDatasets(mode);
      cout << "  -> driver::BuildDatasets returned N=" << datasets.size() << "\n";

      if (datasets.empty())
      {
        cout << ANSI_BOLD_RED << "[FATAL] No datasets selected. Check toggles." << ANSI_RESET << "\n";
        return 1;
      }

      cout << "  Datasets:\n";
      for (const auto& ds : datasets)
      {
        cout << "    - label=" << ds.label
             << "  isSim=" << (ds.isSim ? "true" : "false")
             << "  trigger=" << ds.trigger
             << "  topDir=" << ds.topDirName
             << "  inFile=" << ds.inFilePath
             << "  outBase=" << ds.outBase
             << "\n";
      }

      // ---------------------------------------------------------------------------
      // Open datasets (fail-fast)
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 3] Open all datasets (fail-fast)\n" << ANSI_RESET;

      for (auto& ds : datasets)
      {
        cout << "  -> Opening dataset: " << ANSI_BOLD_YEL << ds.label << ANSI_RESET << "\n";
        cout << "     file=" << ds.inFilePath << "\n";
        cout << "     topDir=" << ds.topDirName << "\n";

        if (!driver::OpenDataset(ds))
        {
          cout << ANSI_BOLD_RED
               << "  [FATAL] Failed to open dataset: " << ds.label << "\n"
               << "         Closing any datasets already opened...\n"
               << ANSI_RESET;

          for (auto& d2 : datasets) driver::CloseDataset(d2);
          return 1;
        }

        cout << ANSI_BOLD_GRN << "     [OK] Opened.\n" << ANSI_RESET;
      }

      // ---------------------------------------------------------------------------
      // Sections 1–3 (per-dataset): event-level QA, failures, isolation
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 4] Sections 1–3 (per-dataset)\n" << ANSI_RESET;

      for (auto& ds : datasets)
      {
        cout << ANSI_BOLD_YEL
             << "\n[DATASET] " << ds.label
             << "  isSim=" << (ds.isSim ? "true" : "false")
             << "  trigger=" << ds.trigger
             << ANSI_RESET << "\n";

        cout << "  Paths:\n"
             << "    inFile   = " << ds.inFilePath << "\n"
             << "    outBase  = " << ds.outBase << "\n"
             << "    topDir   = " << ds.topDirName << "\n";

        cout << "  -> Reading event count...\n";
        const double NevtTotal = ReadEventCount(ds);
        cout << "     NevtTotal (cnt_" << ds.topDirName << " bin1) = "
             << std::fixed << std::setprecision(0) << NevtTotal << "\n";

        cout << "  -> [Section 1] Event-level QA...\n";
        analysis::RunEventLevelQA(ds);
        cout << "     [OK] Event-level QA complete.\n";

        cout << "  -> [triggerAna] doNotScale max-cluster-energy turn-on plots (DATA only)...\n";
        analysis::RunTriggerAna_DoNotScaleMaxClusterEnergy(ds);
        cout << "     [OK] triggerAna complete.\n";

        cout << "  -> [Section 2] Preselection failure table...\n";
        analysis::RunPreselectionFailureTable(ds);
        cout << "     [OK] Preselection failure table complete.\n";

        cout << "  -> [Section 3] General isolation QA...\n";
        analysis::RunIsolationQA(ds);
        cout << "     [OK] Isolation QA complete.\n";
      }

      // ---------------------------------------------------------------------------
      // Section 4: ABCD purity + sideband subtraction with optional leakage factors
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 5] Section 4: ABCD purity + sideband subtraction\n" << ANSI_RESET;
      cout << "  Policy (preserved legacy behavior):\n"
           << "    - Leakage factors are read from SIM if available.\n"
           << "    - Leakage factors are applied ONLY to DATA.\n"
           << "    - SIM runs uncorrected.\n";

      LeakageFactors leakageFromSim;

      cout << "  -> Searching datasets for first SIM source of leakage factors...\n";
      for (auto& ds : datasets)
      {
        if (ds.isSim && !leakageFromSim.available)
        {
          cout << "     Attempt load leakage factors from SIM dataset: " << ds.label << "\n";
          analysis::LoadLeakageFactorsFromSIM(ds, leakageFromSim);
          cout << "     leakageFromSim.available = " << (leakageFromSim.available ? "true" : "false") << "\n";
        }
      }

      if (!leakageFromSim.available)
      {
        cout << ANSI_BOLD_YEL
             << "  [WARN] No leakage factors found from SIM. All datasets will run uncorrected.\n"
             << ANSI_RESET;
      }
      else
      {
        cout << ANSI_BOLD_GRN
             << "  [OK] Leakage factors loaded from SIM. Will apply to DATA only.\n"
             << ANSI_RESET;
      }

      for (auto& ds : datasets)
      {
        LeakageFactors lfForThis;

        const bool applyLeakage = (!ds.isSim && leakageFromSim.available);
        if (applyLeakage)
        {
          lfForThis = leakageFromSim;
          cout << "  -> Running ABCD on DATA dataset: " << ds.label
               << "  (leakage correction: ON)\n";
        }
        else
        {
          lfForThis.available = false;
          cout << "  -> Running ABCD on dataset: " << ds.label
               << "  (leakage correction: OFF)\n";
        }

        analysis::RunABCDPurityAndSidebandSubtraction(ds, lfForThis);
        cout << "     [OK] ABCD purity + sideband subtraction complete.\n";
      }

      // ---------------------------------------------------------------------------
      // Sections 5A–5H: Jet QA suite per dataset (MatchCache per dataset)
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 6] Sections 5A–5H: Jet QA suite (per dataset)\n" << ANSI_RESET;

      map<string, MatchCache> matchCaches;

      for (auto& ds : datasets)
      {
        cout << ANSI_BOLD_YEL << "\n[DATASET JET-QA] " << ds.label << ANSI_RESET << "\n";

        MatchCache& mc = matchCaches[ds.label];
        InitMatchCache(mc);

        cout << "  -> [5A] General jet QA...\n";
        analysis::RunGeneralJetQA(ds);
        cout << "     [OK]\n";

        cout << "  -> [5C] #gamma-jet match QA...\n";
        analysis::RunMatchQA(ds, mc);
        cout << "     [OK]\n";

        cout << "  -> [5D] Selected jet QA...\n";
        analysis::RunSelectedJetQA(ds);
        cout << "     [OK]\n";

        cout << "  -> [5E] xJ + alpha QA...\n";
        analysis::RunXJAlphaQA(ds);
        cout << "     [OK]\n";

        cout << "  -> [5F] JES3 QA...\n";
        analysis::RunJES3QA(ds);
        cout << "     [OK]\n";

        cout << "  -> [5G] Maps QA...\n";
        analysis::RunMapsQA(ds);
        cout << "     [OK]\n";

        cout << "  -> [5H] Unfolding QA...\n";
        analysis::RunUnfoldingQA(ds);
        cout << "     [OK]\n";
      }

      // ---------------------------------------------------------------------------
      // Close + summary
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 7] Summary + close\n" << ANSI_RESET;

      cout << ANSI_BOLD_CYN << "\n[SUMMARY] Histogram coverage (requested vs in-file)" << ANSI_RESET << "\n";

      auto IsHistLikeClass = [](const string& cls) -> bool
      {
        // Inventory class names are like: TH1F, TH2D, TH3F, TProfile, TProfile2D, TProfile3D, ...
        return (cls.rfind("TH", 0) == 0) || (cls.rfind("TProfile", 0) == 0);
      };

      auto PrintTopList = [&](const vector<string>& v, const string& title, int nShow)
      {
        if (v.empty()) return;
        cout << "    " << title << " (showing up to " << nShow << "):\n";
        const int n = std::min((int)v.size(), nShow);
        for (int i = 0; i < n; ++i)
        {
          cout << "      - " << v[i] << "\n";
        }
      };

      for (auto& ds : datasets)
      {
        const string missRawPath = JoinPath(ds.outBase, "missing_hists_" + ds.label + ".txt");
        const string missTabPath = JoinPath(ds.outBase, "missing_hists_" + ds.label + "_TABULATED.txt");
        const string unusedPath  = JoinPath(ds.outBase, "unused_hists_" + ds.label + ".txt");
        const string reqPath     = JoinPath(ds.outBase, "requested_hists_" + ds.label + ".txt");

        // --- Requested (unique fullpaths) ---
        std::set<string> requested;
        for (const auto& kv : ds.requestCounts) requested.insert(kv.first);

        // --- Inventory everything under topDir, then keep only histogram-like objects ---
        vector<InvItem> items;
        CollectInventoryRecursive(ds.topDir, ds.topDirName + "/", items);

        std::set<string> available;
        map<string, InvItem> invByPath;
        for (const auto& it : items)
        {
          if (!IsHistLikeClass(it.cls)) continue;
          available.insert(it.path);
          invByPath[it.path] = it;
        }

        // --- Requested-but-not-in-file (called by code, missing from ROOT file) ---
        vector<string> reqMissing;
        reqMissing.reserve(requested.size());
        for (const auto& p : requested)
        {
          if (available.find(p) == available.end()) reqMissing.push_back(p);
        }

        // --- In-file-but-unused (present in ROOT, never requested by code) ---
        vector<string> unused;
        unused.reserve(available.size());
        for (const auto& p : available)
        {
          if (requested.find(p) == requested.end()) unused.push_back(p);
        }

        // --- Write: requested list (what the code tried to read) ---
        {
          vector<string> lines;
          lines.push_back("path\tcalls");
          for (const auto& p : requested)
          {
            const int calls = (ds.requestCounts.count(p) ? ds.requestCounts.at(p) : 0);
            std::ostringstream s;
            s << p << "\t" << calls;
            lines.push_back(s.str());
          }
          WriteTextFile(reqPath, lines);
        }

        // --- Write: tabulated missing summary ---
        //     Section 1: requested-but-not-in-file  (true missing)
        //     Section 2: logged issues that DO exist in the file (e.g., ZERO_ENTRIES, WRONG_TYPE)
        {
          vector<string> lines;

          lines.push_back("[REQUESTED_BUT_NOT_IN_FILE]");
          lines.push_back("path\tcalls\tloggedCount\tlastReason");
          for (const auto& p : reqMissing)
          {
            const int calls  = (ds.requestCounts.count(p) ? ds.requestCounts.at(p) : 0);
            const int logged = (ds.missingCounts.count(p) ? ds.missingCounts.at(p) : 0);

            string reason = "NOT_IN_FILE";
            auto itR = ds.missingReason.find(p);
            if (itR != ds.missingReason.end()) reason = itR->second;

            std::ostringstream s;
            s << p << "\t" << calls << "\t" << logged << "\t" << reason;
            lines.push_back(s.str());
          }

          lines.push_back("");
          lines.push_back("[LOGGED_ISSUES_PRESENT_IN_FILE]");
          lines.push_back("path\tcalls\tloggedCount\tlastReason");
          for (const auto& kv : ds.missingCounts)
          {
            const string& p = kv.first;

            // If it isn't in-file, it belongs to the section above
            if (available.find(p) == available.end()) continue;

            const int calls  = (ds.requestCounts.count(p) ? ds.requestCounts.at(p) : 0);
            const int logged = kv.second;

            string reason = "(unknown)";
            auto itR = ds.missingReason.find(p);
            if (itR != ds.missingReason.end()) reason = itR->second;

            std::ostringstream s;
            s << p << "\t" << calls << "\t" << logged << "\t" << reason;
            lines.push_back(s.str());
          }

          WriteTextFile(missTabPath, lines);
        }

        // --- Write: unused list (in-file but never requested) ---
        {
          vector<string> lines;
          lines.push_back("path\tclass\tentries");
          for (const auto& p : unused)
          {
            auto it = invByPath.find(p);
            if (it == invByPath.end())
            {
              lines.push_back(p + "\t(n/a)\t(n/a)");
              continue;
            }

            std::ostringstream ent;
            if (it->second.entries < 0) ent << "n/a";
            else ent << std::fixed << std::setprecision(0) << it->second.entries;

            std::ostringstream s;
            s << it->second.path << "\t" << it->second.cls << "\t" << ent.str();
            lines.push_back(s.str());
          }
          WriteTextFile(unusedPath, lines);
        }

        // --- Terminal summary (clear + actionable) ---
        cout << ANSI_BOLD_CYN << "\n  [" << ds.label << "] Histogram I/O coverage" << ANSI_RESET << "\n";
        cout << "    topDir                 : " << ds.topDirName << "\n";
        cout << "    requested (unique)     : " << requested.size()
             << "   (full list: " << reqPath << ")\n";
        cout << "    in file (unique hists) : " << available.size() << "\n";
        cout << "    requested but missing  : " << reqMissing.size()
             << "   (tabulated: " << missTabPath << ")\n";
        cout << "    in file but unused     : " << unused.size()
             << "   (tabulated: " << unusedPath << ")\n";
        cout << "    missing-object log     : occurrences=" << ds.missingCount
             << "  uniqueLogged=" << ds.missingCounts.size()
             << "  rawLog=" << missRawPath << "\n";

        PrintTopList(reqMissing, "Requested-but-not-in-file", 10);
        PrintTopList(unused, "In-file-but-unused", 10);
      }

      cout << "\n  -> Closing datasets...\n";
      for (auto& ds : datasets)
      {
        cout << "     closing: " << ds.label << "\n";
        driver::CloseDataset(ds);
      }

      cout << ANSI_BOLD_CYN << "\nDone.\n" << ANSI_RESET;
      return 0;
    }

} // namespace ARJ


// =============================================================================
// ROOT macro entrypoint
// =============================================================================
int AnalyzeRecoilJets()
{
  return ARJ::Run();
}
