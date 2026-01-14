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
          outPath = JoinPath(ds.outBase, "GeneralEventLevelQA/" + ds.trigger + "/zvtx_DATA_" + ds.trigger + ".png");

        vector<string> lines;
        lines.push_back(TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data());
        DrawAndSaveTH1_Common(ds, hFixed, outPath, "v_{z} [cm]", "Counts", lines, false);

        delete hFixed;
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
          //   DATA -> <outBase>/PurityABCD/<trigger>/preselection/
          // ---------------------------------------------------------------------------
          string outDir;
          if (ds.isSim) outDir = JoinPath(ds.outBase, "PurityABCD/preselection");
          else          outDir = JoinPath(ds.outBase, "PurityABCD/" + ds.trigger + "/preselection");
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
            cTbl.Divide(3,3, 0.002, 0.002);

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
      }

      void RunIsolationQA(Dataset& ds)
      {
        cout << ANSI_BOLD_CYN << "\n==============================\n"
             << "[SECTION 3] GENERAL ISOLATION QA (" << ds.label << ")\n"
             << "==============================" << ANSI_RESET << "\n";

        string outDir;
        if (ds.isSim) outDir = JoinPath(ds.outBase, "isoQAgeneral");
        else          outDir = JoinPath(ds.outBase, "isoQAgeneral/" + ds.trigger);

        EnsureDir(outDir);
        for (const auto& b : PtBins()) EnsureDir(JoinPath(outDir, b.folder));

        // 1) RECO/DATA isolation QA (always runs)
        RunIsolationQA_Reco(ds, outDir);

        // 2) Keep exactly where you had it before (after reco plots, before truth block)
        PrintIsoDecisionTable(ds);

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
        c.Divide(3,3, 0.002, 0.002);

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
          c.Divide(3,3, 0.001, 0.001);

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
        else          outDir = JoinPath(ds.outBase, "PurityABCD/" + ds.trigger);
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

        // 5A: pTcut_noFiducial (all)
        PlotJetQA_AllOrIncl(ds, JoinPath(ds.outBase, "GeneralJetQA/pTcut_noFiducial"),
                           "all", true, false);

        // 5B: fiducial inclusive jets (incl) + event-level
        const string baseIncl = JoinPath(ds.outBase, "GeneralJetQA/pTcutFiducialJets");
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
            else          baseOut = JoinPath(ds.outBase, "RecoilJetQA/MatchQA/" + ds.trigger);
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
                  "c_tbl_dphi", 1500, 1200
                );
                c.Divide(3,3, 0.001, 0.001);

                vector<TH1*> keep;
                keep.reserve(kNPtBins);

                for (int i = 0; i < kNPtBins; ++i)
                {
                  c.cd(i+1);
                  gPad->SetLeftMargin(0.14);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.14);
                  gPad->SetTopMargin(0.10);
                  gPad->SetLogy(logy);

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

              vector<double> x(kNPtBins, 0.0), y02(kNPtBins, 0.0), y04(kNPtBins, 0.0);
              for (int i = 0; i < kNPtBins; ++i)
              {
                const PtBin& b = PtBins()[i];
                x[i] = 0.5*(b.lo + b.hi);

                const double l02 = mc.NphoLead["r02"][i];
                const double l04 = mc.NphoLead["r04"][i];
                const double m02 = mc.NphoMatched["r02"][i];
                const double m04 = mc.NphoMatched["r04"][i];

                y02[i] = (l02 > 0.0) ? (m02/l02) : 0.0;
                y04[i] = (l04 > 0.0) ? (m04/l04) : 0.0;
              }

              // Matched fraction overlay (existing)
              {
                TCanvas c("c_ov_match","c_ov_match",900,700);
                ApplyCanvasMargins1D(c);

                TGraph g02(kNPtBins, &x[0], &y02[0]);
                TGraph g04(kNPtBins, &x[0], &y04[0]);
                g02.SetLineWidth(2); g04.SetLineWidth(2);
                g02.SetMarkerStyle(20); g04.SetMarkerStyle(24);
                g02.SetLineColor(1); g04.SetLineColor(2);
                g02.SetMarkerColor(1); g04.SetMarkerColor(2);

                g02.Draw("ALP");
                g02.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                g02.GetYaxis()->SetTitle("Matched fraction");
                g02.SetMinimum(0.0); g02.SetMaximum(1.05);
                g04.Draw("LP same");

                TLegend leg(0.62,0.78,0.92,0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.033);
                leg.AddEntry(&g02, "r02 (R=0.2)", "lp");
                leg.AddEntry(&g04, "r04 (R=0.4)", "lp");
                leg.Draw();

                DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                DrawLatexLines(0.14,0.78, {"Overlay: Matched fraction vs p_{T}^{#gamma}"}, 0.030, 0.040);

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
            : JoinPath(ds.outBase, "RecoilJetQA/SelectedJetQA/" + ds.trigger);

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
              : JoinPath(ds.outBase, "RecoilJetQA/xJAlpha/" + ds.trigger);

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
              // Only relevant if the user selected the merged photonJet10+20 SIM sample.
              if (CurrentSimSample() != SimSample::kPhotonJet10And20Merged) return true;

              const std::string hay =
                s.label + " " + s.topDirName + " " + s.inFilePath;

              return ContainsAny(hay, {"Merged","merged","10and20","10And20","photon10and20","photonJet10and20","Photon10and20"});
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

            auto RunPair = [&](Dataset& dataDs, Dataset& simDs)
            {
              const std::string dataOutDir =
                dataDs.isSim
                  ? JoinPath(dataDs.outBase, "RecoilJetQA/JES3")
                  : JoinPath(dataDs.outBase, "RecoilJetQA/JES3/" + dataDs.trigger);

              const std::string insBase = JoinPath(dataOutDir, "inSituResidualCalibration");
              EnsureDir(insBase);

              cout << ANSI_BOLD_CYN
                   << "\n============================================================\n"
                   << "[JES3 InSitu] Residual calibration (DATA vs SIM, reco xJ peak)\n"
                   << "  DATA: " << dataDs.label << "  trigger=" << dataDs.trigger << "\n"
                   << "  SIM : " << simDs.label << "  (" << SimSampleLabel(CurrentSimSample()) << ")\n"
                   << "  Output -> " << insBase << "\n"
                   << "============================================================\n"
                   << ANSI_RESET;

              for (const auto& rKey : kRKeys)
              {
                const double R = RFromKey(rKey);

                TH3* hData = GetObj<TH3>(dataDs, "h_JES3_pT_xJ_alpha_" + rKey, false, false, false);
                TH3* hMC   = GetObj<TH3>(simDs,  "h_JES3_pT_xJ_alpha_" + rKey, false, false, false);

                if (!hData || !hMC)
                {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] [JES3 InSitu] Missing reco JES3 xJ TH3 for rKey=" << rKey
                       << "  (DATA=" << (hData ? "OK" : "MISSING")
                       << ", SIM="  << (hMC   ? "OK" : "MISSING") << ")"
                       << ANSI_RESET << "\n";
                  continue;
                }

                const TAxis* axD = hData->GetXaxis();
                const TAxis* axM = hMC->GetXaxis();
                const int nPtD = axD ? axD->GetNbins() : 0;
                const int nPtM = axM ? axM->GetNbins() : 0;
                const int nPt  = std::min(nPtD, nPtM);

                if (nPt <= 0)
                {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] [JES3 InSitu] nPt<=0 for rKey=" << rKey << ANSI_RESET << "\n";
                  continue;
                }
                if (nPtD != nPtM)
                {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] [JES3 InSitu] pT axis mismatch for rKey=" << rKey
                       << "  DATA nPt=" << nPtD << "  SIM nPt=" << nPtM
                       << "  -> using nPt=" << nPt
                       << ANSI_RESET << "\n";
                }

                const std::string rDir   = JoinPath(insBase, rKey);
                const std::string perDir = JoinPath(rDir, "perPtBin");
                EnsureDir(rDir);
                EnsureDir(perDir);

                cout << ANSI_BOLD_CYN
                     << "\n[JES3 InSitu peak-fit table] rKey=" << rKey
                     << " (R=" << std::fixed << std::setprecision(1) << R << ")\n"
                     << ANSI_RESET;

                const int wPt   = 16;
                const int wN    = 12;
                const int wMu   = 10;
                const int wSig  = 10;
                const int wCorr = 12;

                cout << std::left  << std::setw(wPt)   << "pTgamma bin"
                     << std::right << std::setw(wN)    << "N(data)"
                     << std::setw(wMu)   << "muD"
                     << std::setw(wSig)  << "sigD"
                     << std::setw(wN)    << "N(MC)"
                     << std::setw(wMu)   << "muMC"
                     << std::setw(wSig)  << "sigMC"
                     << std::setw(wCorr) << "MC/Data"
                     << "\n";
                cout << std::string(wPt + 2*wN + 2*wMu + 2*wSig + wCorr, '-') << "\n";

                std::vector<std::string> sum;
                sum.push_back(std::string("JES3 in-situ residual (DATA vs SIM)"));
                sum.push_back(std::string("DATA: ") + dataDs.label + "  trigger=" + dataDs.trigger);
                sum.push_back(std::string("SIM : ") + simDs.label);
                sum.push_back(std::string("rKey: ") + rKey + TString::Format("  R=%.1f", R).Data());
                sum.push_back("");
                sum.push_back("Recipe:");
                sum.push_back("  H_i(xJ) = ProjY( h_JES3_pT_xJ_alpha | X=pTgamma bin i, integrate all alpha )");
                sum.push_back("  Peak proxy R_i from iterative Gaussian:");
                sum.push_back("    x0=max-bin center; fit1 [x0±0.15]; refit2 [mu1±1.5*sigma1]; R_i=mu2");
                sum.push_back("  Residual correction factor (apply to DATA jets):  C_i = (R_i^MC / R_i^DATA)");
                sum.push_back("");

                std::vector<double> ptCenters(nPt, 0.0);
                std::vector<double> muD(nPt, 0.0),  muM(nPt, 0.0);
                std::vector<double> sdD(nPt, 0.0),  sdM(nPt, 0.0);
                std::vector<double> cf (nPt, 0.0);

                for (int ib = 1; ib <= nPt; ++ib)
                {
                  const std::string ptLab = AxisBinLabel(axD, ib, "GeV", 0);
                  const double xC = 0.5 * (axD->GetBinLowEdge(ib) + axD->GetBinUpEdge(ib));
                  ptCenters[ib-1] = xC;

                  TH1* hd = ProjectY_AtXbin_TH3(hData, ib,
                    TString::Format("insitu_xJ_data_%s_%d_%s", rKey.c_str(), ib, dataDs.label.c_str()).Data()
                  );
                  TH1* hm = ProjectY_AtXbin_TH3(hMC,   ib,
                    TString::Format("insitu_xJ_mc_%s_%d_%s",   rKey.c_str(), ib, simDs.label.c_str()).Data()
                  );

                  if (hd) { hd->SetDirectory(nullptr); EnsureSumw2(hd); }
                  if (hm) { hm->SetDirectory(nullptr); EnsureSumw2(hm); }

                  const double nData = hd ? hd->GetEntries() : 0.0;
                  const double nMC   = hm ? hm->GetEntries() : 0.0;

                  PeakFitResult fD = FitPeakIterativeGaussian(hd);
                  PeakFitResult fM = FitPeakIterativeGaussian(hm);

                  muD[ib-1] = fD.ok ? fD.mu    : 0.0;
                  sdD[ib-1] = fD.ok ? fD.sigma : 0.0;
                  muM[ib-1] = fM.ok ? fM.mu    : 0.0;
                  sdM[ib-1] = fM.ok ? fM.sigma : 0.0;

                  cf[ib-1]  = (fD.ok && fM.ok) ? SafeDivide(fM.mu, fD.mu, 0.0) : 0.0;

                  cout << std::left  << std::setw(wPt) << ptLab
                       << std::right << std::setw(wN)  << std::fixed << std::setprecision(0) << nData
                       << std::setw(wMu) << std::fixed << std::setprecision(4) << muD[ib-1]
                       << std::setw(wSig)<< std::fixed << std::setprecision(4) << sdD[ib-1]
                       << std::setw(wN)  << std::fixed << std::setprecision(0) << nMC
                       << std::setw(wMu) << std::fixed << std::setprecision(4) << muM[ib-1]
                       << std::setw(wSig)<< std::fixed << std::setprecision(4) << sdM[ib-1]
                       << std::setw(wCorr)<< std::fixed << std::setprecision(4) << cf[ib-1]
                       << "\n";

                  sum.push_back(TString::Format(
                    "pT=%s  Ndata=%.0f muD=%.6f sigD=%.6f  Nmc=%.0f muMC=%.6f sigMC=%.6f  C=muMC/muD=%.6f",
                    ptLab.c_str(),
                    nData, muD[ib-1], sdD[ib-1],
                    nMC,   muM[ib-1], sdM[ib-1],
                    cf[ib-1]
                  ).Data());

                  if (hd)
                  {
                    std::vector<std::string> lines = {
                      "In-situ JES3 residual (DATA): x_{J#gamma}",
                      TString::Format("DATA: %s", dataDs.label.c_str()).Data(),
                      TString::Format("SIM : %s", simDs.label.c_str()).Data(),
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
                      TString::Format("Peak fit: mu=%.4f  sigma=%.4f", muD[ib-1], sdD[ib-1]).Data(),
                      TString::Format("Fit2 range: [%.3f, %.3f]", fD.r2Lo, fD.r2Hi).Data()
                    };

                    DrawXJWithFit(
                      dataDs, hd, fD,
                      JoinPath(perDir, TString::Format("xJ_data_fit_pTbin%d.png", ib).Data()),
                      "x_{J#gamma}", lines,
                      1, 4
                    );
                  }

                  if (hm)
                  {
                    std::vector<std::string> lines = {
                      "In-situ JES3 residual (SIM/MC): x_{J#gamma}",
                      TString::Format("DATA: %s", dataDs.label.c_str()).Data(),
                      TString::Format("SIM : %s", simDs.label.c_str()).Data(),
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
                      TString::Format("Peak fit: mu=%.4f  sigma=%.4f", muM[ib-1], sdM[ib-1]).Data(),
                      TString::Format("Fit2 range: [%.3f, %.3f]", fM.r2Lo, fM.r2Hi).Data()
                    };

                    DrawXJWithFit(
                      simDs, hm, fM,
                      JoinPath(perDir, TString::Format("xJ_mc_fit_pTbin%d.png", ib).Data()),
                      "x_{J#gamma}", lines,
                      2, 4
                    );
                  }

                  if (hd && hm)
                  {
                    std::vector<std::string> lines = {
                      "In-situ JES3 residual overlay (shape): DATA vs MC",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
                      TString::Format("mu(DATA)=%.4f  mu(MC)=%.4f  C=MC/DATA=%.4f", muD[ib-1], muM[ib-1], cf[ib-1]).Data()
                    };

                    DrawOverlayDataVsMCShape(
                      dataDs,
                      hd, hm,
                      JoinPath(perDir, TString::Format("overlay_data_vs_mc_shape_pTbin%d.png", ib).Data()),
                      lines
                    );
                  }

                  if (hd) delete hd;
                  if (hm) delete hm;
                }

                WriteTextFile(JoinPath(rDir, "summary_inSituResidual.txt"), sum);

                // Graphs
                {
                  TCanvas c1(TString::Format("c_mu_vs_pt_%s_%s", dataDs.label.c_str(), rKey.c_str()).Data(), "c_mu_vs_pt", 900, 700);
                  ApplyCanvasMargins1D(c1);

                  TGraph gD(nPt, &ptCenters[0], &muD[0]);
                  TGraph gM(nPt, &ptCenters[0], &muM[0]);

                  gD.SetLineWidth(2);
                  gD.SetMarkerStyle(20);
                  gD.SetLineColor(1);
                  gD.SetMarkerColor(1);

                  gM.SetLineWidth(2);
                  gM.SetMarkerStyle(24);
                  gM.SetLineColor(2);
                  gM.SetMarkerColor(2);

                  gD.Draw("ALP");
                  gD.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                  gD.GetYaxis()->SetTitle("Peak x_{J#gamma} (Gaussian #mu)");
                  gM.Draw("LP same");

                  TLegend leg(0.62, 0.78, 0.92, 0.90);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.033);
                  leg.AddEntry(&gD, "DATA: R_{i} = #mu", "lp");
                  leg.AddEntry(&gM, "MC  : R_{i} = #mu", "lp");
                  leg.Draw();

                  DrawLatexLines(0.14, 0.92, DefaultHeaderLines(dataDs), 0.034, 0.045);
                  DrawLatexLines(0.14, 0.78,
                    {
                      "In-situ JES3 residual calibration",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      TString::Format("DATA=%s   SIM=%s", dataDs.label.c_str(), simDs.label.c_str()).Data()
                    },
                    0.030, 0.040
                  );

                  SaveCanvas(c1, JoinPath(rDir, "graph_peakMu_vs_pTgamma_DATA_vs_MC.png"));
                }

                {
                  TCanvas c2(TString::Format("c_cf_vs_pt_%s_%s", dataDs.label.c_str(), rKey.c_str()).Data(), "c_cf_vs_pt", 900, 700);
                  ApplyCanvasMargins1D(c2);

                  TGraph gC(nPt, &ptCenters[0], &cf[0]);
                  gC.SetLineWidth(2);
                  gC.SetMarkerStyle(20);
                  gC.Draw("ALP");

                  gC.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                  gC.GetYaxis()->SetTitle("Residual correction C_{i} = R_{i}^{MC} / R_{i}^{DATA}");

                  DrawLatexLines(0.14, 0.92, DefaultHeaderLines(dataDs), 0.034, 0.045);
                  DrawLatexLines(0.14, 0.78,
                    {
                      "In-situ JES3 residual correction factor",
                      TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                      "Apply to DATA jet p_{T}:  p_{T}^{jet,corr} = C_{i} #times p_{T}^{jet}"
                    },
                    0.030, 0.040
                  );

                  SaveCanvas(c2, JoinPath(rDir, "graph_residualCorrection_MCOverDATA_vs_pTgamma.png"));
                }

                {
                  std::vector<std::string> corr;
                  corr.push_back("pTgamma_bin  ptCenter  mu_DATA  mu_MC  C=MC/DATA");
                  for (int i = 0; i < nPt; ++i)
                  {
                    const std::string ptLab = AxisBinLabel(axD, i+1, "GeV", 0);
                    corr.push_back(TString::Format(
                      "%s  %.3f  %.6f  %.6f  %.6f",
                      ptLab.c_str(), ptCenters[i], muD[i], muM[i], cf[i]
                    ).Data());
                  }
                  WriteTextFile(JoinPath(rDir, "correctionFactors_MCOverDATA.txt"), corr);
                }
              }
            };

            // Run on every DATA dataset not yet processed, once SIM is available.
            for (Dataset* d : C.data)
            {
              if (!d) continue;
              if (C.doneDataLabels.count(d->label)) continue;

              RunPair(*d, *C.sim);
              C.doneDataLabels.insert(d->label);
            }
        }
  
  
        // =============================================================================
        // JES3 r02 vs r04 overlays (SIM only)
        //
        // Implements the full functionality of the previous RunR02R04Overlays lambda:
        //   - TRUTH xJ integrated over alpha
        //   - RECO  xJ integrated over alpha
        //   - RECO  xJ with alpha cuts
        //
        // Output under:
        //   <outDir>/r02_r04/...
        //
        // Colors and styling preserved.
        // =============================================================================
        void JES3_R02R04Overlays_MaybeRun(Dataset& ds, const std::string& outDir)
        {
            if (!ds.isSim) return;

            const std::string ovBase = JoinPath(outDir, "r02_r04");
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
              [&](const TH3* h02, const TH3* h04,
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
              const int nPt = std::min(n02, n04);

              // Per-bin overlay PNGs
              for (int ib = 1; ib <= nPt; ++ib)
              {
                TH1* a = nullptr;
                TH1* b = nullptr;

                if (useAlphaCut)
                {
                  a = ProjectY_AtXbin_AndAlphaMax_TH3(h02, ib, alphaMax,
                        TString::Format("xJ_ov_r02_alphaLT%.2f_b%d", alphaMax, ib).Data());
                  b = ProjectY_AtXbin_AndAlphaMax_TH3(h04, ib, alphaMax,
                        TString::Format("xJ_ov_r04_alphaLT%.2f_b%d", alphaMax, ib).Data());
                }
                else
                {
                  a = ProjectY_AtXbin_TH3(h02, ib, TString::Format("xJ_ov_r02_int_b%d", ib).Data());
                  b = ProjectY_AtXbin_TH3(h04, ib, TString::Format("xJ_ov_r04_int_b%d", ib).Data());
                }

                if (!a && !b) { if (a) delete a; if (b) delete b; continue; }

                const std::string ptLab = AxisBinLabel(h02->GetXaxis(), ib, "GeV", 0);

                TCanvas c(TString::Format("c_ov_%s_%d", outDirHere.c_str(), ib).Data(), "c_ov", 900, 700);
                ApplyCanvasMargins1D(c);

                // Style: force clean, filled markers and matching error-bar colors
                if (a)
                {
                  a->SetDirectory(nullptr);
                  EnsureSumw2(a);

                  a->SetTitle("");
                  a->SetLineWidth(2);
                  a->SetLineColor(2);      // red
                  a->SetMarkerStyle(20);   // filled circle
                  a->SetMarkerSize(1.05);
                  a->SetMarkerColor(2);
                  a->SetFillStyle(0);

                    a->GetXaxis()->SetTitle(xTitle.c_str());
                    a->GetYaxis()->SetTitle((ds.isSim && bothPhoton10and20sim) ? "Counts / pb^{-1}" : "Counts");
                }
                if (b)
                {
                  b->SetDirectory(nullptr);
                  EnsureSumw2(b);

                  b->SetTitle("");
                  b->SetLineWidth(2);
                  b->SetLineColor(4);      // blue
                  b->SetMarkerStyle(20);   // filled circle
                  b->SetMarkerSize(1.05);
                  b->SetMarkerColor(4);
                  b->SetFillStyle(0);

                    b->GetXaxis()->SetTitle(xTitle.c_str());
                    b->GetYaxis()->SetTitle((ds.isSim && bothPhoton10and20sim) ? "Counts / pb^{-1}" : "Counts");
                }

                TH1* first  = a ? a : b;
                TH1* second = (first == a) ? b : a;

                double ymax = 0.0;
                if (a) ymax = std::max(ymax, a->GetMaximum());
                if (b) ymax = std::max(ymax, b->GetMaximum());
                if (first) first->SetMaximum(ymax * 1.25);

                // Draw with error bars + markers
                if (first)  first->Draw("E1");
                if (second) second->Draw("E1 same");

                // Legend should reflect markers+errors
                TLegend leg(0.70, 0.78, 0.92, 0.90);
                leg.SetTextFont(42);
                leg.SetTextSize(0.035);
                if (a) leg.AddEntry(a, "r02 (red)", "ep");
                if (b) leg.AddEntry(b, "r04 (blue)", "ep");
                leg.Draw();

                // Header stays at top-left
                DrawLatexLines(0.14, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);

                // Info block just under the header
                std::vector<std::string> lines = headerLines;
                lines.push_back(TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());
                if (useAlphaCut) lines.push_back(TString::Format("#alpha < %.2f", alphaMax).Data());
                DrawLatexLines(0.14, 0.875, lines, 0.030, 0.038);

                SaveCanvas(c, JoinPath(outDirHere, TString::Format("overlay_pTbin%d.png", ib).Data()));

                if (a) delete a;
                if (b) delete b;
              }

              // 3x3 table(s) of overlays (linear y only)
              const int perPage = 9;
              int page = 0;

              for (int start = 1; start <= nPt; start += perPage)
              {
                ++page;

                TCanvas c(
                  TString::Format("c_tbl_%s_p%d", outDirHere.c_str(), page).Data(),
                  "c_tbl_overlay", 1500, 1200
                );
                c.Divide(3,3, 0.001, 0.001);

                std::vector<TH1*> keep;
                keep.reserve(2 * perPage);

                for (int k = 0; k < perPage; ++k)
                {
                  const int ib = start + k;
                  c.cd(k+1);

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

                  TH1* a = nullptr;
                  TH1* b = nullptr;

                  if (useAlphaCut)
                  {
                    a = ProjectY_AtXbin_AndAlphaMax_TH3(h02, ib, alphaMax,
                          TString::Format("tbl_r02_alphaLT%.2f_p%d_b%d", alphaMax, page, ib).Data());
                    b = ProjectY_AtXbin_AndAlphaMax_TH3(h04, ib, alphaMax,
                          TString::Format("tbl_r04_alphaLT%.2f_p%d_b%d", alphaMax, page, ib).Data());
                  }
                  else
                  {
                    a = ProjectY_AtXbin_TH3(h02, ib, TString::Format("tbl_r02_int_p%d_b%d", page, ib).Data());
                    b = ProjectY_AtXbin_TH3(h04, ib, TString::Format("tbl_r04_int_p%d_b%d", page, ib).Data());
                  }

                  if (!a && !b)
                  {
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextSize(0.06);
                    t.DrawLatex(0.15, 0.55, "MISSING");
                    continue;
                  }

                  // Style (must set marker color or ROOT will look black in small pads)
                  if (a)
                  {
                    a->SetDirectory(nullptr);
                    EnsureSumw2(a);

                    a->SetTitle("");
                    a->SetLineWidth(2);
                    a->SetLineColor(2);
                    a->SetMarkerStyle(20);
                    a->SetMarkerSize(0.95);
                    a->SetMarkerColor(2);
                    a->SetFillStyle(0);

                      a->GetXaxis()->SetTitle(xTitle.c_str());
                      a->GetYaxis()->SetTitle((ds.isSim && bothPhoton10and20sim) ? "Counts / pb^{-1}" : "Counts");
                  }
                  if (b)
                  {
                    b->SetDirectory(nullptr);
                    EnsureSumw2(b);

                    b->SetTitle("");
                    b->SetLineWidth(2);
                    b->SetLineColor(4);
                    b->SetMarkerStyle(20);
                    b->SetMarkerSize(0.95);
                    b->SetMarkerColor(4);
                    b->SetFillStyle(0);

                      b->GetXaxis()->SetTitle(xTitle.c_str());
                      b->GetYaxis()->SetTitle((ds.isSim && bothPhoton10and20sim) ? "Counts / pb^{-1}" : "Counts");
                  }

                  TH1* first  = a ? a : b;
                  TH1* second = (first == a) ? b : a;

                  double ymax = 0.0;
                  if (a) ymax = std::max(ymax, a->GetMaximum());
                  if (b) ymax = std::max(ymax, b->GetMaximum());
                  if (first) first->SetMaximum(ymax * 1.25);

                  if (first)  first->Draw("E1");
                  if (second) second->Draw("E1 same");

                  // Big centered title per pad (kept as in your helper)
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


                  // Legend top-right (requested); no overlay text
                  TLegend leg(0.8, 0.75, 0.88, 0.9);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.055);
                  leg.SetFillStyle(0);
                  leg.SetBorderSize(0);
                  if (a) leg.AddEntry(a, "r02", "ep");
                  if (b) leg.AddEntry(b, "r04", "ep");
                  leg.DrawClone();

                  if (useAlphaCut)
                  {
                    TLatex talpha;
                    talpha.SetNDC(true);
                    talpha.SetTextFont(42);
                    talpha.SetTextAlign(13);
                    talpha.SetTextSize(0.045);
                    talpha.DrawLatex(0.16, 0.86, TString::Format("#alpha < %.2f", alphaMax).Data());
                  }

                  if (a) keep.push_back(a);
                  if (b) keep.push_back(b);
                }

                std::string outName;
                if (nPt <= perPage)
                {
                  outName = useAlphaCut
                    ? TString::Format("table3x3_overlay_alphaLT%.2f.png", alphaMax).Data()
                    : "table3x3_overlay_integratedAlpha.png";
                }
                else
                {
                  outName = useAlphaCut
                    ? TString::Format("table3x3_overlay_alphaLT%.2f_page%d.png", alphaMax, page).Data()
                    : TString::Format("table3x3_overlay_integratedAlpha_page%d.png", page).Data();
                }

                SaveCanvas(c, JoinPath(outDirHere, outName));

                for (auto* h : keep) delete h;
              }
            };

            // -------------------- TRUTH (reco-conditioned): integrated alpha --------------------
            TH3* hTr02 = GetObj<TH3>(ds, "h_JES3Truth_pT_xJ_alpha_r02", true, true, true);
            TH3* hTr04 = GetObj<TH3>(ds, "h_JES3Truth_pT_xJ_alpha_r04", true, true, true);

            if (hTr02 && hTr04)
            {
              const std::string ovTruthRecoCond = JoinPath(ovBase, "xJ_truth_integratedAlpha_recoConditioned");
              DrawOverlayPair_TH3xJ(
                hTr02, hTr04,
                ovTruthRecoCond,
                "x_{J#gamma}^{truth}",
                {"TRUTH (reco-conditioned): x_{J#gamma}^{truth}", "Overlay: r02 red, r04 blue"},
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
              const std::string ovTruthPure = JoinPath(ovBase, "xJ_truth_integratedAlpha_pureTruth");
              DrawOverlayPair_TH3xJ(
                hTrPure02, hTrPure04,
                ovTruthPure,
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

            // -------------------- RECO (existing): integrated alpha --------------------
            TH3* hRe02 = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_r02", true, true, true);
            TH3* hRe04 = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_r04", true, true, true);

            if (hRe02 && hRe04)
            {
              const std::string ovReco = JoinPath(ovBase, "xJ_reco_integratedAlpha");
              DrawOverlayPair_TH3xJ(
                hRe02, hRe04,
                ovReco,
                "x_{J#gamma}",
                {"RECO: x_{J#gamma}", "Overlay: r02 red, r04 blue"},
                false, 0.0
              );

              // -------------------- RECO (existing): alpha-cut overlays --------------------
              const std::vector<double> alphaMaxCuts = {0.20, 0.30, 0.40, 0.50};
              const std::string ovRecoCutsBase = JoinPath(ovBase, "xJ_reco_alphaCuts");
              EnsureDir(ovRecoCutsBase);

              for (double aMax : alphaMaxCuts)
              {
                const std::string aDir = JoinPath(ovRecoCutsBase, AlphaTag(aMax));
                DrawOverlayPair_TH3xJ(
                  hRe02, hRe04,
                  aDir,
                  "x_{J#gamma}",
                  {"RECO: x_{J#gamma}", "Overlay: r02 red, r04 blue"},
                  true, aMax
                );
              }
            }
            else
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] RECO overlays skipped: missing h_JES3_pT_xJ_alpha_r02 or r04 in dataset " << ds.label
                   << ANSI_RESET << "\n";
            }
        }

        void RunJES3QA(Dataset& ds)
        {
            cout << ANSI_BOLD_CYN << "\n==============================\n"
                 << "[SECTION 5F] JES3 QA (3D maps) (" << ds.label << ")\n"
                 << "==============================" << ANSI_RESET << "\n";

            string outDir = ds.isSim
              ? JoinPath(ds.outBase, "RecoilJetQA/JES3")
              : JoinPath(ds.outBase, "RecoilJetQA/JES3/" + ds.trigger);

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
              [&](TH3* hReco_xJ, TH3* hReco_j1, TH3* hTrut_xJ, TH3* hTrut_j1, TH3* hTrutPure_xJ)->const TAxis*
            {
              return
                (hReco_xJ ? hReco_xJ->GetXaxis()
                : (hReco_j1 ? hReco_j1->GetXaxis()
                : (hTrut_xJ ? hTrut_xJ->GetXaxis()
                : (hTrut_j1 ? hTrut_j1->GetXaxis()
                : (hTrutPure_xJ ? hTrutPure_xJ->GetXaxis() : nullptr)))));
            };

            struct Jes3Dirs
            {
              string rOut;
              string dir2D;
              string dirSumm;
              string dirXJProj;
              string dirXJProjReco;
              string dirXJProjTruthRecoCond;
              string dirXJProjTruthPure;
            };

            auto MakeJes3Dirs =
              [&](const string& rKey, double /*R*/)->Jes3Dirs
            {
              Jes3Dirs D;
              D.rOut     = JoinPath(outDir, rKey);
              D.dir2D    = JoinPath(D.rOut, "2DMaps");
              D.dirSumm  = JoinPath(D.rOut, "summaries");

              // xJ(=Y) projections after integrating over alpha(=Z)
              D.dirXJProj               = JoinPath(D.rOut, "xJ_fromJES3");
              D.dirXJProjReco           = JoinPath(D.dirXJProj, "RECO");
              D.dirXJProjTruthRecoCond  = JoinPath(D.dirXJProj, "TRUTH_recoConditioned");
              D.dirXJProjTruthPure      = JoinPath(D.dirXJProj, "TRUTH_pure");

              EnsureDir(D.rOut);
              EnsureDir(D.dir2D);
              EnsureDir(D.dirSumm);
              EnsureDir(D.dirXJProj);
              EnsureDir(D.dirXJProjReco);
              EnsureDir(D.dirXJProjTruthRecoCond);
              EnsureDir(D.dirXJProjTruthPure);

              return D;
            };

            struct Jes3Hists
            {
              TH3* hReco_xJ      = nullptr;
              TH3* hTrut_xJ      = nullptr;  // reco-conditioned truth (existing)
              TH3* hTrutPure_xJ  = nullptr;  // pure-truth xJgamma TH3 (no reco gating)
              TH3* hReco_j1      = nullptr;
              TH3* hTrut_j1      = nullptr;
            };

            auto LoadJes3Hists =
              [&](const string& rKey)->Jes3Hists
            {
              Jes3Hists H;
              H.hReco_xJ      = GetObj<TH3>(ds, "h_JES3_pT_xJ_alpha_" + rKey, true, true, true);
              H.hTrut_xJ      = GetObj<TH3>(ds, "h_JES3Truth_pT_xJ_alpha_" + rKey, true, true, true);
              H.hTrutPure_xJ  = GetObj<TH3>(ds, "h_JES3TruthPure_pT_xJ_alpha_" + rKey, true, true, true);
              H.hReco_j1      = GetObj<TH3>(ds, "h_JES3_pT_jet1Pt_alpha_" + rKey, true, true, true);
              H.hTrut_j1      = GetObj<TH3>(ds, "h_JES3Truth_pT_jet1Pt_alpha_" + rKey, true, true, true);
              return H;
            };

            auto HasAnyJes3 =
              [&](const Jes3Hists& H)->bool
            {
              return (H.hReco_xJ || H.hTrut_xJ || H.hTrutPure_xJ || H.hReco_j1 || H.hTrut_j1);
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
              const TAxis* axPt = PickPtAxis(H.hReco_xJ, H.hReco_j1, H.hTrut_xJ, H.hTrut_j1, H.hTrutPure_xJ);
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
                vector<string> lines = {
                  "JES3 (RECO): x_{J#gamma}",
                  TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                  TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data(),
                  "Projection: integrated over #alpha",
                  "Filled with: (p_{T}^{#gamma}, x_{J#gamma}, #alpha)"
                };

                DrawAndSaveTH1_Common(ds, xJ_re,
                  JoinPath(D.dirXJProjReco, TString::Format("xJ_reco_integratedAlpha_pTbin%d.png", ib).Data()),
                  "x_{J#gamma}", "Counts", lines, false, false, 0.0, "E1");

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
                    "JES3 (TRUTH reco-conditioned): x_{J#gamma}^{truth}",
                    TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                    TString::Format("p_{T}^{#gamma,truth}: %s", ptLab.c_str()).Data(),
                    "Filled with: (tPt, xJt=tj1Pt/tPt, aT=tj2Pt/tPt)"
                  };

                  DrawAndSaveTH1_Common(ds, xJ_tr,
                    JoinPath(D.dirXJProjTruthRecoCond, TString::Format("xJ_truth_integratedAlpha_recoConditioned_pTbin%d.png", ib).Data()),
                    "x_{J#gamma}^{truth}", "Counts", lines, false, false, 0.0, "E1");
                };

                auto SaveXJTruthPurePNGs =
                  [&](int ib, const string& ptLab)
                {
                  if (!H.hTrutPure_xJ) return;

                  TH1* xJ_pure = ProjectY_AtXbin_TH3(
                    H.hTrutPure_xJ, ib,
                    TString::Format("jes3_xJ_trPure_%s_%d", rKey.c_str(), ib).Data()
                  );

                  if (!xJ_pure) return;

                  xJ_pure->SetName(TString::Format("h_xJ_trPure_%s_pTbin%d", rKey.c_str(), ib).Data());

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

                  // TRUTH reco-conditioned (existing h_JES3Truth_pT_xJ_alpha_<rKey>)
                  SaveXJTruthPNGs(P.xJ_tr, ib, P.ptLab);

                  // TRUTH pure (new h_JES3TruthPure_pT_xJ_alpha_<rKey>)
                  SaveXJTruthPurePNGs(ib, P.ptLab);

                // RECO/TRUTH xJ-alpha maps and jet1Pt-alpha maps
                SaveJes3Maps2D_ForBin(ib, P.ptLab);

                CleanupBinPack(P);
              }

              // -------------------------------------------------------------------------
              // 3x3 tables + text summary (existing behavior preserved)
              // -------------------------------------------------------------------------

              auto Make3x3Table_xJ_FromTH3 =
                [&](const TH3* h3, const string& outBaseDir, const string& tag, bool logy)
              {
                if (!h3) return;

                const int n = h3->GetXaxis()->GetNbins();
                const int perPage = 9;

                int page = 0;
                for (int start = 1; start <= n; start += perPage)
                {
                  ++page;

                  TCanvas c(
                    TString::Format("c_tbl_xJ_%s_%s_%s_%s_p%d",
                      ds.label.c_str(),
                      rKey.c_str(),
                      tag.c_str(),
                      logy ? "logy" : "lin",
                      page).Data(),
                    "c_tbl_xJ", 1500, 1200
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
                    gPad->SetLogy(logy);

                    if (ib > n)
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
                      hx->GetYaxis()->SetTitle((ds.isSim && bothPhoton10and20sim) ? "Counts / pb^{-1}" : "Counts");

                    if (logy)
                    {
                      const double minPos = SmallestPositiveBinContent(hx);
                      hx->SetMinimum((minPos > 0.0) ? (0.5 * minPos) : 1e-6);
                    }

                    hx->Draw("E1");

                    const string ptLab = AxisBinLabel(h3->GetXaxis(), ib, "GeV", 0);

                    vector<string> lines;
                    lines.push_back(TString::Format("JES3 %s: x_{J#gamma} ( #int d#alpha )", tag.c_str()).Data());
                    lines.push_back(TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data());
                    lines.push_back(TString::Format("p_{T}^{#gamma}: %s", ptLab.c_str()).Data());
                    if (logy) lines.push_back("Y-axis: log scale");

                    DrawLatexLines(0.16, 0.90, lines, 0.040, 0.050);

                    keep.push_back(hx);
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
                      hx->GetYaxis()->SetTitle((ds.isSim && bothPhoton10and20sim) ? "Counts / pb^{-1}" : "Counts");
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
              Make3x3Table_xJ_FromTH3(H.hReco_xJ,     D.dirXJProjReco,          "RECO",  false);
              Make3x3Table_xJ_FromTH3(H.hTrut_xJ,     D.dirXJProjTruthRecoCond, "TRUTH", false);
              Make3x3Table_xJ_FromTH3(H.hTrutPure_xJ, D.dirXJProjTruthPure,     "TRUTH", false);

              Make3x3Table_xJ_FromTH3(H.hReco_xJ,     D.dirXJProjReco,          "RECO",  true);
              Make3x3Table_xJ_FromTH3(H.hTrut_xJ,     D.dirXJProjTruthRecoCond, "TRUTH", true);
              Make3x3Table_xJ_FromTH3(H.hTrutPure_xJ, D.dirXJProjTruthPure,     "TRUTH", true);

              // RECO alpha-cut tables into subfolders (requested)
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

                for (double aMax : alphaMaxCuts)
                {
                  const string aDir = JoinPath(dirAlphaBase, AlphaTag(aMax));
                  EnsureDir(aDir);
                  Make3x3Table_xJ_FromTH3_AlphaCut(H.hReco_xJ, aDir, "RECO", aMax);
                }
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

            JES3_R02R04Overlays_MaybeRun(ds, outDir);
            
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
            : JoinPath(ds.outBase, "RecoilJetQA/Maps/" + ds.trigger);
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
              : JoinPath(ds.outBase, "RecoilJetQA/Unfolding/" + ds.trigger);

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

            auto LoadUnfoldHists = [&](const string& rKey)->UnfoldHists
            {
              UnfoldHists H;

              H.hReco   = GetObj<TH2>(ds, "h2_unfoldReco_pTgamma_xJ_incl_" + rKey, true, true, true);
              H.hTruth  = GetObj<TH2>(ds, "h2_unfoldTruth_pTgamma_xJ_incl_" + rKey, true, true, true);
              H.hFakes  = GetObj<TH2>(ds, "h2_unfoldRecoFakes_pTgamma_xJ_incl_" + rKey, true, true, true);
              H.hMisses = GetObj<TH2>(ds, "h2_unfoldTruthMisses_pTgamma_xJ_incl_" + rKey, true, true, true);
              H.hResp   = GetObj<TH2>(ds, "h2_unfoldResponse_pTgamma_xJ_incl_" + rKey, true, true, true);

              H.hRecoDphi  = GetObj<TH2>(ds, "h2_unfoldReco_pTgamma_dphi_incl_" + rKey, true, true, true);
              H.hTruthDphi = GetObj<TH2>(ds, "h2_unfoldTruth_pTgamma_dphi_incl_" + rKey, true, true, true);

              return H;
            };

            auto HasAnyUnfold = [&](const UnfoldHists& H)->bool
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
              Save2D(rOut, H.hReco,   "unfold_reco_pTgamma_vs_xJ.png",
                     "p_{T}^{#gamma,reco} [GeV]",  "x_{J#gamma}^{reco}",  "Counts",
                     {string("Unfold input: RECO counts"), rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              Save2D(rOut, H.hTruth,  "unfold_truth_pTgamma_vs_xJ.png",
                     "p_{T}^{#gamma,truth} [GeV]", "x_{J#gamma}^{truth}", "Counts",
                     {string("Unfold truth: TRUTH counts"), rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              Save2D(rOut, H.hFakes,  "unfold_recoFakes_pTgamma_vs_xJ.png",
                     "p_{T}^{#gamma,reco} [GeV]",  "x_{J#gamma}^{reco}",  "Counts (fakes)",
                     {string("Unfold input: RECO fakes"), rKey + TString::Format(" (R=%.1f)", R).Data()},
                     true);

              // Response matrix (SCAT) (kept exactly as your block)
              if (H.hResp)
              {
                TH2* hc = CloneTH2(H.hResp, "resp_scatter_clone");
                if (hc)
                {
                  const int oldOptStat = gStyle->GetOptStat();
                  gStyle->SetOptStat(1110);
                  hc->SetStats(1);

                  TCanvas c("c_resp_scatter","c_resp_scatter",950,780);
                  ApplyCanvasMargins2D(c);

                  hc->SetTitle("");
                  hc->GetXaxis()->SetTitle("global bin (truth: p_{T}^{#gamma}, x_{J})");
                  hc->GetYaxis()->SetTitle("global bin (reco: p_{T}^{#gamma}, x_{J})");
                  hc->GetZaxis()->SetTitle("Counts");

                  hc->SetMarkerStyle(6);
                  hc->SetMarkerSize(0.6);
                  hc->SetMarkerColor(kBlue+1);

                  hc->Draw("SCAT");

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
              EnsureDir(dphiTruthDir);
              EnsureDir(dphiProjDir);

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
                const int perPage = 9;
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

                const int nPtUse = std::min((int)PtBins().size(), H.hRecoDphi->GetXaxis()->GetNbins());
                for (int i = 0; i < nPtUse; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  const int xbin = i + 1;

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
                    "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1");

                  TH1* ps = CloneTH1(p, TString::Format("p_absDphi_reco_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data());
                  NormalizeVisible(ps);

                  DrawAndSaveTH1_Common(ds, ps,
                    JoinPath(perDir, "absDphi_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", lines, false, false, 0.0, "E1");

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

                const int nPtUse = std::min((int)PtBins().size(), H.hTruthDphi->GetXaxis()->GetNbins());
                for (int i = 0; i < nPtUse; ++i)
                {
                  const PtBin& b = PtBins()[i];
                  const int xbin = i + 1;

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
                    "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1");

                  TH1* ps = CloneTH1(p, TString::Format("p_absDphi_truth_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), i).Data());
                  NormalizeVisible(ps);

                  DrawAndSaveTH1_Common(ds, ps,
                    JoinPath(perDir, "absDphi_shape.png"),
                    "|#Delta#phi| [rad]", "A.U.", lines, false, false, 0.0, "E1");

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
              const int nPt = std::min(nPtReco, nPtTruth);

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

              vector<double> xCenters(nPt, 0.0);
              vector<double> effVsPt(nPt, 0.0);
              vector<double> purVsPt(nPt, 0.0);

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

              // Per pT bin: terminal line + text summary + overlay truth vs reco xJ shapes
              for (int ib = 1; ib <= nPt; ++ib)
              {
                const double xC = 0.5 * (H.hTruth->GetXaxis()->GetBinLowEdge(ib) + H.hTruth->GetXaxis()->GetBinUpEdge(ib));
                xCenters[ib-1] = xC;

                const string ptLabTruth = AxisBinLabel(H.hTruth->GetXaxis(), ib, "GeV", 0);
                const string ptLabReco  = AxisBinLabel(H.hReco ->GetXaxis(),  ib, "GeV", 0);

                const double Ntruth = H.hTruth->Integral(ib, ib, 0, H.hTruth->GetYaxis()->GetNbins() + 1);
                const double Nreco  = H.hReco ->Integral(ib, ib, 0, H.hReco ->GetYaxis()->GetNbins() + 1);

                const double Nmiss  = (H.hMisses ? H.hMisses->Integral(ib, ib, 0, H.hMisses->GetYaxis()->GetNbins() + 1) : 0.0);
                const double Nfake  = (H.hFakes  ? H.hFakes ->Integral(ib, ib, 0, H.hFakes ->GetYaxis()->GetNbins() + 1) : 0.0);

                const double eff = SafeDivide(Ntruth - Nmiss, Ntruth, 0.0);
                const double pur = SafeDivide(Nreco  - Nfake, Nreco,  0.0);

                effVsPt[ib-1] = eff;
                purVsPt[ib-1] = pur;

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

                // Existing truth-vs-reco xJ shape overlay per pT bin (kept)
                {
                  TH1D* pxTruth = H.hTruth->ProjectionY(TString::Format("pxTruth_%s_%d", rKey.c_str(), ib).Data(), ib, ib);
                  TH1D* pxReco  = H.hReco ->ProjectionY(TString::Format("pxReco_%s_%d", rKey.c_str(), ib).Data(),  ib, ib);
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

                    TCanvas c(TString::Format("c_unf_ov_%s_%d", rKey.c_str(), ib).Data(), "c_unf_ov", 900,700);
                    ApplyCanvasMargins1D(c);

                    const double maxv = std::max(pxTruth->GetMaximum(), pxReco->GetMaximum());
                    pxTruth->SetMaximum(maxv * 1.25);
                    pxTruth->SetTitle("");
                    pxTruth->GetXaxis()->SetTitle("x_{J}");
                    pxTruth->GetYaxis()->SetTitle("A.U.");
                    pxTruth->Draw("hist");
                    pxReco->Draw("hist same");

                    TLegend leg(0.60,0.78,0.92,0.90);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.033);
                    leg.AddEntry(pxTruth, "Truth (shape)", "l");
                    leg.AddEntry(pxReco,  "Reco (shape)",  "l");
                    leg.Draw();

                    DrawLatexLines(0.14,0.92, DefaultHeaderLines(ds), 0.034, 0.045);
                    DrawLatexLines(0.14,0.78, {string("Truth vs reco x_{J} shape"), rKey, TString::Format("p_{T}^{#gamma}: %s", ptLabTruth.c_str()).Data()}, 0.030, 0.040);

                    SaveCanvas(c, JoinPath(rOut, TString::Format("overlay_truth_vs_reco_xJ_shape_pTbin%d.png", ib).Data()));

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
                TCanvas c1(TString::Format("c_eff_vs_pt_%s", rKey.c_str()).Data(), "c_eff_vs_pt", 900,700);
                ApplyCanvasMargins1D(c1);
                TGraph g(nPt, &xCenters[0], &effVsPt[0]);
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
              {
                TCanvas c2(TString::Format("c_pur_vs_pt_%s", rKey.c_str()).Data(), "c_pur_vs_pt", 900,700);
                ApplyCanvasMargins1D(c2);
                TGraph g(nPt, &xCenters[0], &purVsPt[0]);
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
    inline bool OpenDataset(Dataset& ds)
    {
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
          ds.outBase    = kOutPPBase;
          datasets.push_back(std::move(ds));
        }

        return datasets;
    }

    inline bool MaybeBuildMergedSIM(RunMode mode)
    {
        // Only relevant when a SIM-including mode is running AND the merged SIM sample was selected.
        if (mode == RunMode::kPPDataOnly) return true;

        const SimSample ss = CurrentSimSample();
        if (ss != SimSample::kPhotonJet10And20Merged) return true;

        const bool ok = BuildMergedSIMFile_Photon10And20(
          kInSIM10,
          kInSIM20,
          kMergedSIMOut,
          kDirSIM,
          6692.7611,
          105.79868
        );

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
           << "    isPPdataOnly         = " << (isPPdataOnly ? "true" : "false") << "\n"
           << "    isSimAndDataPP       = " << (isSimAndDataPP ? "true" : "false") << "\n"
           << "    isPhotonJet10        = " << (isPhotonJet10 ? "true" : "false") << "\n"
           << "    isPhotonJet20        = " << (isPhotonJet20 ? "true" : "false") << "\n"
           << "    bothPhoton10and20sim = " << (bothPhoton10and20sim ? "true" : "false") << "\n";

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
      // Optional SIM slice merge
      // ---------------------------------------------------------------------------
      cout << ANSI_BOLD_CYN << "\n[STEP 1] Optional SIM slice merge (photonJet10 + photonJet20)\n" << ANSI_RESET;
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
