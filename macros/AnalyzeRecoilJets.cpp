// AnalyzeRecoilJets.cpp
//
// Usage (ROOT):
//   root -l -q AnalyzeRecoilJets.cpp
// -----------------------------------------------------------------------------

#include "AnalyzeRecoilJets.h"
#include "TGraphAsymmErrors.h"

#ifdef __has_include
#if __has_include("RooUnfoldResponse.h")
#define ARJ_HAVE_ROOUNFOLD 1
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfold.h"
#else
#define ARJ_HAVE_ROOUNFOLD 0
#endif
#else
#define ARJ_HAVE_ROOUNFOLD 0
#endif

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
        if (!ds.centSuffix.empty())
        {
          vector<InvItem> filtered;
          filtered.reserve(items.size());
          for (const auto& it : items)
          {
            const bool isCentHist = (it.path.find("_cent_") != string::npos);
            if (isCentHist && it.path.find(ds.centSuffix) == string::npos) continue;
            filtered.push_back(it);
          }
          items.swap(filtered);
        }
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
          if (thr == 10) return kViolet+1;
          if (thr == 12) return kCyan+2;
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

        std::vector<LoadedPair> pairs;

        TIter next(ds.file->GetListOfKeys());
        while (TKey* key = dynamic_cast<TKey*>(next()))
        {
          const std::string dirKey = key->GetName();
          if (dirKey.empty()) continue;
          if (dirKey.find(baselineTag) == 0) continue;

          TH1* hProbe = getHist(dirKey);
          if (!hProbe) continue;

          const std::string baselineKey = baselineTag + dirKey;
          TH1* hBase = getHist(baselineKey);
          if (!hBase)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Missing doNotScale baseline histogram for probe " << dirKey << "\n"
                 << "       Need: " << baselineKey << "/" << prefix << baselineKey << "\n"
                 << ANSI_RESET << "\n";
            continue;
          }

          const std::string basisKey = deduceBasisKey(dirKey);
          if (basisKey.empty())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Could not deduce common baseline group for probe " << dirKey
                 << ANSI_RESET << "\n";
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
          pairs.push_back(P);
        }

        if (pairs.empty())
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] No pairwise doNotScale numerator/baseline histogram pairs found in " << ds.inFilePath
               << ANSI_RESET << "\n";
          return;
        }

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

        auto drawGroupTurnOnOverlay = [&](const std::vector<LoadedPair>& loaded, const std::string& groupOutDir, const std::string& filenameTag = "")
        {
          std::vector<TH1*> ratioHists;
          std::vector<double> x95s;
          ratioHists.reserve(loaded.size());
          x95s.reserve(loaded.size());

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
            hRatio->SetMarkerStyle(markerForProbe(P.probeKey));
            hRatio->SetMarkerSize(1.0);
            hRatio->SetMarkerColor(colorForProbe(P.probeKey));
            hRatio->SetLineColor(colorForProbe(P.probeKey));
            hRatio->SetLineWidth(3);

            ratioHists.push_back(hRatio);
            x95s.push_back(findXAtEff(hRatio, 0.95));

            delete hBase;
            delete hProbe;
          }

          if (ratioHists.empty()) return;

          TCanvas c(TString::Format("c_trigAna_groupTurnOn_%s", loaded.front().groupFolder.c_str()).Data(),
                    TString::Format("c_trigAna_groupTurnOn_%s", loaded.front().groupFolder.c_str()).Data(),
                    900, 700);
          c.cd();
          c.SetLeftMargin(0.14);
          c.SetRightMargin(0.05);
          c.SetBottomMargin(0.14);
          c.SetTopMargin(0.08);
          c.SetTicks(1,1);

          const int dsThr = extractPhotonThresholdGeV(ds.trigger);
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

          TLegend leg(0.35, 0.3, 0.79, 0.45);
          leg.SetBorderSize(0);
          leg.SetFillStyle(0);
          leg.SetTextSize(0.024);
          for (std::size_t i = 0; i < loaded.size() && i < ratioHists.size(); ++i)
          {
            if (x95s[i] > 0.0)
            {
              leg.AddEntry(ratioHists[i],
                TString::Format("%s (95%%=%.2f GeV)",
                  loaded[i].probeLabel.c_str(),
                  x95s[i]).Data(),
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

          {
            TLatex tDS;
            tDS.SetNDC(true);
            tDS.SetTextFont(42);
            tDS.SetTextAlign(33);
            tDS.SetTextSize(0.045);
            tDS.DrawLatex(0.93, 0.97, isAuAuOnly ? "AuAu" : (isRun25pp ? "Run25pp" : "Run24pp"));
          }

          const std::string outPng = JoinPath(groupOutDir, "hMaxClusterEnergy_groupTurnOnOverlay" + filenameTag + ".png");
          SaveCanvas(c, outPng);
          cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";

            for (auto* h : ratioHists) delete h;
          };

          auto drawGroupMaxClusterOverlay = [&](const std::vector<LoadedPair>& loaded, const std::string& groupOutDir, const std::string& filenameTag = "")
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
                std::string baseLbl = loaded.front().baselineLabel;
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

            {
              TLatex tDS;
              tDS.SetNDC(true);
              tDS.SetTextFont(42);
              tDS.SetTextAlign(33);
              tDS.SetTextSize(0.045);
              tDS.DrawLatex(0.93, 0.97, isAuAuOnly ? "AuAu" : (isRun25pp ? "Run25pp" : "Run24pp"));
            }

            const std::string outPng = JoinPath(groupOutDir, "hMaxClusterEnergy_groupOverlay" + filenameTag + ".png");
            SaveCanvas(c, outPng);
            cout << ANSI_BOLD_GRN << "[WROTE] " << outPng << ANSI_RESET << "\n";

            delete hBaseGroup;
            for (auto* h : probeHists) delete h;
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

          if (isAuAuOnly && folder == "commonBasis_MBD_NS_geq_2_vtx_lt_150")
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

            if (isAuAuOnly && kv.first == "commonBasis_MBD_NS_geq_2_vtx_lt_150")
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

      void RunPi0QA(Dataset& ds)
      {
        const bool isSimMBDataset = (ds.isSim && (isSimMB || isSimJet5));
        if (ds.isSim && !isSimMBDataset) return;
        if (!isSimMBDataset && ds.trigger != kTriggerPP) return;
        if (!ds.file) return;

        cout << ANSI_BOLD_CYN << "\n==============================\n"
                << "[pi0 QA] corrected vs no-asinh-correction (" << ds.label << ")\n"
                << "==============================" << ANSI_RESET << "\n";

        const string outDir = isSimMBDataset
             ? JoinPath(ds.outBase, "pi0_QA")
             : JoinPath(OutputPP(), "pi0_QA");
            EnsureDir(outDir);

        // SIM: pi0 histograms live under the "SIM" topDir alongside all other SIM hists.
        // DATA: pi0 histograms live under "MBD_NandS_geq_1".
        const string pi0DirName = isSimMBDataset ? "SIM" : "MBD_NandS_geq_1";
        TDirectory* dir = ds.file->GetDirectory(pi0DirName.c_str());
        if (!dir)
        {
          cout << ANSI_BOLD_YEL
                 << "[WARN] Missing pi0 directory: " << pi0DirName
                 << ANSI_RESET << "\n";
          return;
        }

        TH2* h2Corr = dynamic_cast<TH2*>(dir->Get("h2_pi0_mass_vs_leadcluspt_corr"));
        TH2* h2NoCorr = dynamic_cast<TH2*>(dir->Get("h2_pi0_mass_vs_leadcluspt_nocorr"));
        TH2* h2Pi0PtCorr = dynamic_cast<TH2*>(dir->Get("h2_pi0_mass_vs_pi0pt_corr"));
        TH2* h2Pi0PtNoCorr = dynamic_cast<TH2*>(dir->Get("h2_pi0_mass_vs_pi0pt_nocorr"));

        if (!h2Corr || !h2NoCorr)
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] Missing pi0 TH2(s) needed for pi0 QA.\n"
               << "       Need:\n"
               << "         " << pi0DirName << "/h2_pi0_mass_vs_leadcluspt_corr\n"
               << "         " << pi0DirName << "/h2_pi0_mass_vs_leadcluspt_nocorr\n"
               << ANSI_RESET << "\n";
          return;
        }

        const string datasetTitle = isSimMBDataset ? (isSimJet5 ? "InclusiveJet5 SIM" : "MinBias SIM (DETROIT)") : "Run24pp";
        const bool vzIsInteger = (std::fabs(vzCutCm - std::round(vzCutCm)) < 1e-6);
        const string vertexLabel =
          vzIsInteger
            ? string(TString::Format("|v_{z}| < %.0f cm", vzCutCm).Data())
            : string(TString::Format("|v_{z}| < %.1f cm", vzCutCm).Data());

        const int nLeadPtBins = 4;
        const double ptLo[nLeadPtBins] = {1.0, 2.0, 3.0, 5.0};
        const double ptHi[nLeadPtBins] = {2.0, 3.0, 5.0, -1.0};
        const char* ptLabels[nLeadPtBins] = {"1-2", "2-3", "3-5", "> 5"};

        struct Pi0FitResult
        {
          bool ok = false;
          double mean = 0.0;
          double meanErr = 0.0;
          double sigma = 0.0;
          double sigmaErr = 0.0;
          TF1* func = nullptr;
        };

        auto LeadPhotonPtTitle = [&](int i)->string
        {
          if (ptHi[i] > ptLo[i])
          {
            return TString::Format("p_{T}^{#gamma} = %s GeV", ptLabels[i]).Data();
          }
          return "p_{T}^{#gamma} > 5 GeV";
        };

        auto Pi0PtTitle = [&](int i)->string
        {
          if (ptHi[i] > ptLo[i])
          {
            return TString::Format("p_{T}^{#pi^{0}} = %s GeV", ptLabels[i]).Data();
          }
          return "p_{T}^{#pi^{0}} > 5 GeV";
        };

        auto DrawCenteredHeader = [&](const string& qualifier, double size)
        {
          TLatex title;
          title.SetNDC(true);
          title.SetTextFont(42);
          title.SetTextAlign(23);
          title.SetTextSize(size);
          title.DrawLatex(
            0.50, 0.965,
            TString::Format("%s, #pi^{0} with/without b-correction (%s)",
                            datasetTitle.c_str(), qualifier.c_str()).Data()
          );
        };

        auto DrawTopRightSelection = [&](double x, double y, double size, double dy)
        {
          TLatex t;
          t.SetNDC(true);
          t.SetTextFont(42);
          t.SetTextAlign(33);
          t.SetTextSize(size);
          t.DrawLatex(x, y, "Trigger: MBD NS #geq 1");
          t.DrawLatex(x, y - 1.0 * dy, vertexLabel.c_str());
          t.DrawLatex(x, y - 2.0 * dy, "#alpha #leq 0.6");
          t.DrawLatex(x, y - 3.0 * dy, "E_{#gamma} #geq 1 GeV, #chi^{2} #leq 4");
        };

        auto ProjectMass = [&](TH2* h2, double lo, double hi, const string& newName)->TH1D*
        {
          if (!h2) return nullptr;

          const int ny = h2->GetYaxis()->GetNbins();
          int ybinLo = h2->GetYaxis()->FindBin(lo + 1e-6);
          int ybinHi = (hi > lo) ? h2->GetYaxis()->FindBin(hi - 1e-6) : ny;

          if (ybinLo < 1) ybinLo = 1;
          if (ybinHi > ny) ybinHi = ny;
          if (ybinHi < ybinLo) return nullptr;

          TH1D* h = h2->ProjectionX(newName.c_str(), ybinLo, ybinHi);
          if (h) h->SetDirectory(nullptr);
          return h;
        };

        auto FitPi0 = [&](TH1* h, const string& tag)->Pi0FitResult
        {
          Pi0FitResult R;
          if (!h) return R;
          if (h->Integral(0, h->GetNbinsX() + 1) <= 0.0) return R;

          const int bLo = h->GetXaxis()->FindBin(0.08);
          const int bHi = h->GetXaxis()->FindBin(0.20);

          int bMax = h->GetMaximumBin();
          double cMax = -1.0;
          for (int ib = bLo; ib <= bHi; ++ib)
          {
            const double c = h->GetBinContent(ib);
            if (c > cMax)
            {
              cMax = c;
              bMax = ib;
            }
          }

          if (cMax <= 0.0) cMax = h->GetMaximum();
          if (!(cMax > 0.0)) return R;

          const double m0 = h->GetXaxis()->GetBinCenter(bMax);

          TF1 f(TString::Format("f_pi0_%s", tag.c_str()).Data(), "gaus(0)+pol3(3)", 0.05, 0.25);
          f.SetLineWidth(2);
          f.SetNpx(500);
          f.SetParameters(cMax, m0, 0.010, 0.0, 0.0, 0.0, 0.0);
          f.SetParLimits(1, 0.10, 0.17);
          f.SetParLimits(2, 0.003, 0.040);

          h->Fit(&f, "RQ0N");

          const double mean  = f.GetParameter(1);
          const double sigma = std::fabs(f.GetParameter(2));

          if (!std::isfinite(mean) || !std::isfinite(sigma) || !(sigma > 0.0))
          {
            return R;
          }

          R.ok       = true;
          R.mean     = mean;
          R.meanErr  = f.GetParError(1);
          R.sigma    = sigma;
          R.sigmaErr = f.GetParError(2);
          R.func     = static_cast<TF1*>(f.Clone(TString::Format("f_pi0_%s_clone", tag.c_str()).Data()));
          return R;
        };

        auto DrawPi0Table = [&](TH2* h2CorrIn,
                                TH2* h2NoCorrIn,
                                const string& canvasName,
                                const string& outName,
                                const string& tagPrefix,
                                auto&& TitleFunc)
        {
          if (!h2CorrIn || !h2NoCorrIn) return;

          TCanvas cTbl(canvasName.c_str(), canvasName.c_str(), 1800, 1400);
          cTbl.Divide(2, 2, 0.001, 0.001);

          vector<TH1*> keepAlive;
          keepAlive.reserve(2 * nLeadPtBins);

          vector<TF1*> keepFits;
          keepFits.reserve(2 * nLeadPtBins);

          vector<TLegend*> keepLegends;
          keepLegends.reserve(nLeadPtBins);

          double meanCorrLocal[nLeadPtBins];
          double meanNoCorrLocal[nLeadPtBins];
          double resCorrLocal[nLeadPtBins];
          double resNoCorrLocal[nLeadPtBins];
          bool okMeanCorrLocal[nLeadPtBins];
          bool okMeanNoCorrLocal[nLeadPtBins];

          for (int i = 0; i < nLeadPtBins; ++i)
          {
            meanCorrLocal[i] = 0.0;
            meanNoCorrLocal[i] = 0.0;
            resCorrLocal[i] = 0.0;
            resNoCorrLocal[i] = 0.0;
            okMeanCorrLocal[i] = false;
            okMeanNoCorrLocal[i] = false;
          }

          for (int i = 0; i < nLeadPtBins; ++i)
          {
            cTbl.cd(i + 1);

            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.10);
            gPad->SetTicks(1,1);

            TH1D* hCorr = ProjectMass(
              h2CorrIn,
              ptLo[i],
              ptHi[i],
              TString::Format("h_pi0_corr_%s_%d", tagPrefix.c_str(), i).Data()
            );

            TH1D* hNoCorr = ProjectMass(
              h2NoCorrIn,
              ptLo[i],
              ptHi[i],
              TString::Format("h_pi0_nocorr_%s_%d", tagPrefix.c_str(), i).Data()
            );

            const string ptTitle = TitleFunc(i);

            const double iCorr = hCorr ? hCorr->Integral(0, hCorr->GetNbinsX() + 1) : 0.0;
            const double iNoCorr = hNoCorr ? hNoCorr->Integral(0, hNoCorr->GetNbinsX() + 1) : 0.0;

            if (iCorr <= 0.0 && iNoCorr <= 0.0)
            {
              DrawCenteredHeader(ptTitle, 0.036);
              DrawTopRightSelection(0.94, 0.72, 0.040, 0.065);

              TLatex t;
              t.SetNDC(true);
              t.SetTextFont(42);
              t.SetTextAlign(22);
              t.SetTextSize(0.075);
              t.DrawLatex(0.50, 0.55, "MISSING");

              if (hCorr) delete hCorr;
              if (hNoCorr) delete hNoCorr;
              continue;
            }

            if (hCorr)
            {
              hCorr->SetTitle("");
              hCorr->SetLineWidth(2);
              hCorr->SetLineColor(kRed + 1);
              hCorr->SetMarkerStyle(20);
              hCorr->SetMarkerSize(0.85);
              hCorr->SetMarkerColor(kRed + 1);
              hCorr->SetFillStyle(0);
              hCorr->SetMinimum(0.0);
              hCorr->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
              hCorr->GetYaxis()->SetTitle("Counts");
              hCorr->GetXaxis()->SetRangeUser(0.02, 0.30);
              hCorr->GetXaxis()->SetTitleSize(0.055);
              hCorr->GetYaxis()->SetTitleSize(0.055);
              hCorr->GetXaxis()->SetLabelSize(0.045);
              hCorr->GetYaxis()->SetLabelSize(0.045);
              hCorr->GetYaxis()->SetTitleOffset(1.15);
            }

            if (hNoCorr)
            {
              hNoCorr->SetTitle("");
              hNoCorr->SetLineWidth(2);
              hNoCorr->SetLineColor(kBlack);
              hNoCorr->SetMarkerStyle(24);
              hNoCorr->SetMarkerSize(0.85);
              hNoCorr->SetMarkerColor(kBlack);
              hNoCorr->SetFillStyle(0);
              hNoCorr->SetMinimum(0.0);
              hNoCorr->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
              hNoCorr->GetYaxis()->SetTitle("Counts");
              hNoCorr->GetXaxis()->SetRangeUser(0.02, 0.30);
              hNoCorr->GetXaxis()->SetTitleSize(0.055);
              hNoCorr->GetYaxis()->SetTitleSize(0.055);
              hNoCorr->GetXaxis()->SetLabelSize(0.045);
              hNoCorr->GetYaxis()->SetLabelSize(0.045);
              hNoCorr->GetYaxis()->SetTitleOffset(1.15);
            }

            Pi0FitResult fitCorrRes;
            Pi0FitResult fitNoCorrRes;

            if (hCorr)
            {
              fitCorrRes = FitPi0(hCorr, TString::Format("%s_corr_%d", tagPrefix.c_str(), i).Data());
              if (fitCorrRes.ok)
              {
                meanCorrLocal[i] = fitCorrRes.mean;
                resCorrLocal[i] = fitCorrRes.sigma / fitCorrRes.mean;
                okMeanCorrLocal[i] = true;
              }
            }

            if (hNoCorr)
            {
              fitNoCorrRes = FitPi0(hNoCorr, TString::Format("%s_nocorr_%d", tagPrefix.c_str(), i).Data());
              if (fitNoCorrRes.ok)
              {
                meanNoCorrLocal[i] = fitNoCorrRes.mean;
                resNoCorrLocal[i] = fitNoCorrRes.sigma / fitNoCorrRes.mean;
                okMeanNoCorrLocal[i] = true;
              }
            }

            double ymax = 0.0;
            if (hCorr) ymax = std::max(ymax, hCorr->GetMaximum());
            if (hNoCorr) ymax = std::max(ymax, hNoCorr->GetMaximum());

            TH1* first = hCorr ? hCorr : hNoCorr;
            if (first)
            {
              first->SetMinimum(0.0);
              first->SetMaximum((ymax > 0.0) ? (1.35 * ymax) : 1.0);
              first->Draw("E1");
            }

            if (hCorr && hCorr != first) hCorr->Draw("E1 SAME");
            if (hNoCorr && hNoCorr != first) hNoCorr->Draw("E1 SAME");

            if (fitCorrRes.func)
            {
              fitCorrRes.func->SetLineColor(kRed + 1);
              fitCorrRes.func->SetLineWidth(2);
              fitCorrRes.func->SetNpx(500);
              fitCorrRes.func->Draw("SAME");
              keepFits.push_back(fitCorrRes.func);
            }

            if (fitNoCorrRes.func)
            {
              fitNoCorrRes.func->SetLineColor(kBlack);
              fitNoCorrRes.func->SetLineWidth(2);
              fitNoCorrRes.func->SetNpx(500);
              fitNoCorrRes.func->Draw("SAME");
              keepFits.push_back(fitNoCorrRes.func);
            }

            auto* leg = new TLegend(0.16, 0.72, 0.53, 0.82);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.037);
            if (hCorr)   leg->AddEntry(hCorr,   "with b = 0.15", "lep");
            if (hNoCorr) leg->AddEntry(hNoCorr, "no asinh correction", "lep");
            leg->Draw();
            keepLegends.push_back(leg);

            DrawCenteredHeader(ptTitle, 0.036);
            DrawTopRightSelection(0.94, 0.72, 0.040, 0.065);

            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(13);
            t.SetTextSize(0.038);
            if (okMeanCorrLocal[i] && okMeanNoCorrLocal[i] && resNoCorrLocal[i] > 0.0)
            {
              const double resolutionGainPct = 100.0 * (resNoCorrLocal[i] - resCorrLocal[i]) / resNoCorrLocal[i];
              t.SetTextColor(kBlack);
              t.DrawLatex(0.16, 0.87,
                TString::Format("Resolution gain with correction: %.2f%%", resolutionGainPct).Data()
              );
            }
            t.SetTextColor(kBlack);

            if (hCorr) keepAlive.push_back(hCorr);
            if (hNoCorr) keepAlive.push_back(hNoCorr);
          }

          SaveCanvas(cTbl, JoinPath(outDir, outName));

          for (auto* h : keepAlive) delete h;
          for (auto* f : keepFits) delete f;
          for (auto* leg : keepLegends) delete leg;
        };

        double x[nLeadPtBins];
        double ex[nLeadPtBins];
        double meanCorr[nLeadPtBins];
        double meanCorrErr[nLeadPtBins];
        double sigmaCorr[nLeadPtBins];
        double sigmaCorrErr[nLeadPtBins];
        double resCorr[nLeadPtBins];
        double resCorrErr[nLeadPtBins];
        double meanNoCorr[nLeadPtBins];
        double meanNoCorrErr[nLeadPtBins];
        double sigmaNoCorr[nLeadPtBins];
        double sigmaNoCorrErr[nLeadPtBins];
        double resNoCorr[nLeadPtBins];
        double resNoCorrErr[nLeadPtBins];
        bool okMeanCorr[nLeadPtBins];
        bool okMeanNoCorr[nLeadPtBins];

        for (int i = 0; i < nLeadPtBins; ++i)
        {
          x[i] = i + 1.0;
          ex[i] = 0.0;

          meanCorr[i] = 0.0;
          meanCorrErr[i] = 0.0;
          sigmaCorr[i] = 0.0;
          sigmaCorrErr[i] = 0.0;
          resCorr[i] = 0.0;
          resCorrErr[i] = 0.0;

          meanNoCorr[i] = 0.0;
          meanNoCorrErr[i] = 0.0;
          sigmaNoCorr[i] = 0.0;
          sigmaNoCorrErr[i] = 0.0;
          resNoCorr[i] = 0.0;
          resNoCorrErr[i] = 0.0;

          okMeanCorr[i] = false;
          okMeanNoCorr[i] = false;
        }

        DrawPi0Table(
          h2Corr,
          h2NoCorr,
          "c_pi0_table",
          "table2x2_pi0_mass_leadPhotonPt_corr_vs_nocorr.pdf",
          "leadPt",
          LeadPhotonPtTitle
        );

        auto DrawSummaryGraph =
          [&](const string& outName,
              const string& yTitle,
              const string& summaryLine,
              const double yCorrIn[],
              const double yCorrErrIn[],
              const bool okCorrIn[],
              const double yNoCorrIn[],
              const double yNoCorrErrIn[],
              const bool okNoCorrIn[],
              bool forceZeroMin)
        {
          vector<double> xCorr;
          vector<double> exCorr;
          vector<double> yCorr;
          vector<double> eyCorr;

          vector<double> xNoCorr;
          vector<double> exNoCorr;
          vector<double> yNoCorr;
          vector<double> eyNoCorr;

          double yMin = std::numeric_limits<double>::max();
          double yMax = -std::numeric_limits<double>::max();

          for (int i = 0; i < nLeadPtBins; ++i)
          {
            if (okCorrIn[i] && std::isfinite(yCorrIn[i]))
            {
              xCorr.push_back(x[i]);
              exCorr.push_back(0.0);
              yCorr.push_back(yCorrIn[i]);
              eyCorr.push_back(yCorrErrIn[i]);
              yMin = std::min(yMin, yCorrIn[i]);
              yMax = std::max(yMax, yCorrIn[i]);
            }

            if (okNoCorrIn[i] && std::isfinite(yNoCorrIn[i]))
            {
              xNoCorr.push_back(x[i]);
              exNoCorr.push_back(0.0);
              yNoCorr.push_back(yNoCorrIn[i]);
              eyNoCorr.push_back(yNoCorrErrIn[i]);
              yMin = std::min(yMin, yNoCorrIn[i]);
              yMax = std::max(yMax, yNoCorrIn[i]);
            }
          }

          if (!(yMax > -std::numeric_limits<double>::max()))
          {
            yMin = 0.0;
            yMax = 1.0;
          }

          double pad = yMax - yMin;
          if (!(pad > 0.0)) pad = (yMax != 0.0 ? 0.25 * std::fabs(yMax) : 0.25);

          double frameMin = forceZeroMin ? 0.0 : (yMin - 0.20 * pad);
          double frameMax = yMax + 0.25 * pad;

          if (forceZeroMin && frameMax <= 0.0) frameMax = 1.0;
          if (!forceZeroMin && frameMax <= frameMin) frameMax = frameMin + 1.0;

          TCanvas c(TString::Format("c_%s", outName.c_str()).Data(), "c_pi0_summary", 900, 700);
          ApplyCanvasMargins1D(c);
          c.SetTopMargin(0.10);
          c.SetTicks(1,1);

          TH1F hFrame(TString::Format("hFrame_%s", outName.c_str()).Data(), "", nLeadPtBins, 0.5, nLeadPtBins + 0.5);
          hFrame.SetDirectory(nullptr);
          hFrame.SetStats(0);
          hFrame.SetMinimum(frameMin);
          hFrame.SetMaximum(frameMax);
          hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
          hFrame.GetYaxis()->SetTitle(yTitle.c_str());
          hFrame.GetXaxis()->SetTitleSize(0.050);
          hFrame.GetYaxis()->SetTitleSize(0.050);
          hFrame.GetXaxis()->SetLabelSize(0.040);
          hFrame.GetYaxis()->SetLabelSize(0.040);
          hFrame.GetYaxis()->SetTitleOffset(1.20);
          hFrame.GetXaxis()->CenterLabels(true);

          for (int i = 1; i <= nLeadPtBins; ++i)
          {
            hFrame.GetXaxis()->SetBinLabel(i, ptLabels[i - 1]);
          }

          hFrame.Draw();

          TGraphErrors* gCorr = nullptr;
          TGraphErrors* gNoCorr = nullptr;

          if (!xCorr.empty())
          {
            gCorr = new TGraphErrors((int)xCorr.size(), &xCorr[0], &yCorr[0], &exCorr[0], &eyCorr[0]);
            gCorr->SetLineWidth(2);
            gCorr->SetLineColor(kRed + 1);
            gCorr->SetMarkerStyle(20);
            gCorr->SetMarkerSize(1.0);
            gCorr->SetMarkerColor(kRed + 1);
            gCorr->Draw("PE1 SAME");
          }

          if (!xNoCorr.empty())
          {
            gNoCorr = new TGraphErrors((int)xNoCorr.size(), &xNoCorr[0], &yNoCorr[0], &exNoCorr[0], &eyNoCorr[0]);
            gNoCorr->SetLineWidth(2);
            gNoCorr->SetLineColor(kBlack);
            gNoCorr->SetMarkerStyle(24);
            gNoCorr->SetMarkerSize(1.0);
            gNoCorr->SetMarkerColor(kBlack);
            gNoCorr->Draw("PE1 SAME");
          }

          TLegend leg(0.16, 0.74, 0.45, 0.86);
          leg.SetBorderSize(0);
          leg.SetFillStyle(0);
          leg.SetTextSize(0.034);
          if (gCorr)   leg.AddEntry(gCorr,   "with b = 0.15", "lep");
          if (gNoCorr) leg.AddEntry(gNoCorr, "no asinh correction", "lep");
          leg.Draw();

          DrawCenteredHeader(summaryLine, 0.036);
          DrawTopRightSelection(0.93, 0.74, 0.036, 0.060);

          SaveCanvas(c, JoinPath(outDir, outName));

          if (gCorr) delete gCorr;
          if (gNoCorr) delete gNoCorr;
        };

        DrawSummaryGraph(
          "pi0_mean_vs_leadPhotonPt.png",
          "Gaussian mean [GeV]",
          "Gaussian mean vs p_{T}^{#gamma}",
          meanCorr,
          meanCorrErr,
          okMeanCorr,
          meanNoCorr,
          meanNoCorrErr,
          okMeanNoCorr,
          false
        );

        DrawSummaryGraph(
          "pi0_sigma_vs_leadPhotonPt.png",
          "Gaussian #sigma [GeV]",
          "Gaussian #sigma vs p_{T}^{#gamma}",
          sigmaCorr,
          sigmaCorrErr,
          okMeanCorr,
          sigmaNoCorr,
          sigmaNoCorrErr,
          okMeanNoCorr,
          true
        );

        DrawSummaryGraph(
          "pi0_resolution_vs_leadPhotonPt.png",
          "#sigma / #mu",
          "Resolution (#sigma / #mu) vs p_{T}^{#gamma}",
          resCorr,
          resCorrErr,
          okMeanCorr,
          resNoCorr,
          resNoCorrErr,
          okMeanNoCorr,
          true
        );

        if (h2Pi0PtCorr && h2Pi0PtNoCorr)
        {
          DrawPi0Table(
            h2Pi0PtCorr,
            h2Pi0PtNoCorr,
            "c_pi0_table_pi0pt",
            "table2x2_pi0_mass_pi0Pt_corr_vs_nocorr.pdf",
            "pi0Pt",
            Pi0PtTitle
          );
        }
        else
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] Missing pi0 TH2(s) needed for additional pi0-pT QA table.\n"
               << "       Need:\n"
               << "         MBD_NandS_geq_1/h2_pi0_mass_vs_pi0pt_corr\n"
               << "         MBD_NandS_geq_1/h2_pi0_mass_vs_pi0pt_nocorr\n"
               << ANSI_RESET << "\n";
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
               << std::setw(wN) << "weta>=0.6"
               << std::setw(wN) << "et1<=0.6"
               << std::setw(wN) << "et1>=1.0"
               << std::setw(wN) << "et1Out"
               << std::setw(wN) << "e11>=0.98"
               << std::setw(wN) << "e32<=0.8"
               << std::setw(wN) << "e32>=1.0"
               << std::setw(wN) << "e32Out"
               << "\n";
          cout << string(wBin + 8*wN, '-') << "\n";

          // Requested x-axis labels (keep these exact)
          const char* xLabels[8] = {
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H"
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
                gPad->SetLeftMargin(0.16);
                gPad->SetRightMargin(0.04);
                gPad->SetTopMargin(0.09);
                gPad->SetBottomMargin(0.18);
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
            hAxis->GetXaxis()->LabelsOption("h");

            hAxis->GetXaxis()->SetLabelSize(compact ? 0.060 : 0.055);
            hAxis->GetXaxis()->SetLabelOffset(compact ? 0.010 : 0.012);

            // Y-axis title: bigger + closer in the 3x3 table
            hAxis->GetYaxis()->SetTitleSize(compact ? 0.058 : 0.055);
            hAxis->GetYaxis()->SetTitleOffset(compact ? 1.15 : 1.05);
            hAxis->GetYaxis()->SetLabelSize(compact ? 0.046 : 0.045);

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
            t.SetTextSize(compact ? 0.050 : 0.034);


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
              TLatex tt;
              tt.SetTextFont(42);
              tt.SetNDC();
              tt.SetTextAlign(23);
              tt.SetTextSize(0.060);
              tt.DrawLatex(0.50, 0.965,
                TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data()
              );
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
              const char* xLabelsTbl[8] = {
                "A",
                "B",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H"
              };

              const int binColorsTbl[8] = {
                kGray+1,
                kGreen+2,
                kRed+1,
                kMagenta+1,
                kOrange+7,
                kYellow+2,
                kAzure+1,
                kCyan+1
              };

              const int nCols = 3;
              const int nRows = 3;

              TCanvas cTbl(
                  TString::Format("c_preFail_table3x3_%s", ds.label.c_str()).Data(),
                  "c_preFail_table3x3", 3000, 2200
              );
              cTbl.Divide(nCols, nRows, 0.0, 0.0);

              vector<TObject*> keepTbl;
              keepTbl.reserve(kNPtBins * 12);

              const double leftOuter   = 0.036;
              const double rightOuter  = 0.018;
              const double bottomOuter = 0.025;
              const double headerBand  = 0.078;
              const double gapX        = 0.014;
              const double gapY        = 0.026;

              const double padW = (1.0 - leftOuter - rightOuter - (nCols - 1)*gapX) / nCols;
              const double padH = (1.0 - bottomOuter - headerBand - (nRows - 1)*gapY) / nRows;

              for (int i = 0; i < kNPtBins; ++i)
              {
                TPad* pad = dynamic_cast<TPad*>(cTbl.cd(i+1));
                if (!pad) continue;

                const int row = i / nCols;
                const int col = i % nCols;

                const double x1 = leftOuter + col * (padW + gapX);
                const double x2 = x1 + padW;
                const double y2 = 1.0 - headerBand - row * (padH + gapY);
                const double y1 = y2 - padH;

                pad->SetPad(x1, y1, x2, y2);
                pad->SetLeftMargin(0.13);
                pad->SetRightMargin(0.030);
                pad->SetTopMargin(0.085);
                pad->SetBottomMargin(0.140);
                pad->SetTicks(1,1);

                double vv[8];
                for (int ib = 0; ib < 8; ++ib) vv[ib] = valsByPt[i][ib];

                double ymax = 0.0;
                for (int ib = 0; ib < 8; ++ib) ymax = std::max(ymax, vv[ib]);
                const double yMaxPlot = (ymax > 0.0) ? (1.30 * ymax) : 1.0;

                TH1F* hAxis = new TH1F(
                  TString::Format("h_preFailTblAxis_%s_%s", ds.label.c_str(), bins[i].folder.c_str()).Data(),
                  "",
                  8, 0.5, 8.5
                );
                hAxis->SetDirectory(nullptr);
                hAxis->SetStats(0);
                hAxis->SetMinimum(0.0);
                hAxis->SetMaximum(yMaxPlot);

                for (int ib = 1; ib <= 8; ++ib) hAxis->GetXaxis()->SetBinLabel(ib, xLabelsTbl[ib-1]);

                hAxis->GetYaxis()->SetTitle("Fail counts");
                hAxis->GetXaxis()->SetTitle("");
                hAxis->GetXaxis()->LabelsOption("h");
                hAxis->GetXaxis()->SetLabelFont(62);
                hAxis->GetXaxis()->SetLabelSize(0.060);
                hAxis->GetXaxis()->SetLabelOffset(0.008);

                hAxis->GetYaxis()->SetTitleSize(0.048);
                hAxis->GetYaxis()->SetTitleOffset(1.15);
                hAxis->GetYaxis()->SetLabelSize(0.040);

                hAxis->SetLineColor(1);
                hAxis->SetLineWidth(2);
                hAxis->SetFillStyle(0);
                hAxis->Draw("hist");

                for (int ib = 1; ib <= 8; ++ib)
                {
                  TH1F* hb = new TH1F(
                    TString::Format("h_preFailTblBar_%s_%s_b%d", ds.label.c_str(), bins[i].folder.c_str(), ib).Data(),
                    "",
                    8, 0.5, 8.5
                  );
                  hb->SetDirectory(nullptr);
                  hb->SetStats(0);

                  hb->SetBinContent(ib, vv[ib-1]);
                  hb->SetFillStyle(1001);
                  hb->SetFillColor(binColorsTbl[ib-1]);
                  hb->SetLineColor(1);
                  hb->SetLineWidth(2);

                  hb->SetBarWidth(0.88);
                  hb->SetBarOffset(0.06);

                  hb->Draw("BAR SAME");
                  keepTbl.push_back(hb);
                }

                TLatex t;
                t.SetTextFont(42);
                t.SetTextAlign(22);
                t.SetTextSize(0.043);

                for (int ib = 1; ib <= 8; ++ib)
                {
                  const double y = vv[ib-1];
                  if (y <= 0.0) continue;

                  const double x = hAxis->GetXaxis()->GetBinCenter(ib);
                  const double yText = std::min(y + 0.022*yMaxPlot, 0.93*yMaxPlot);
                  t.DrawLatex(x, yText, TString::Format("%.0f", y).Data());
                }

                TLatex tt;
                tt.SetTextFont(42);
                tt.SetNDC();
                tt.SetTextAlign(23);
                tt.SetTextSize(0.052);
                tt.DrawLatex(0.50, 0.958,
                  TString::Format("p_{T}^{#gamma}: %d-%d GeV", bins[i].lo, bins[i].hi).Data()
                );

                keepTbl.push_back(hAxis);
              }

              cTbl.cd(0);

              TLatex tGlobal;
              tGlobal.SetNDC();
              tGlobal.SetTextFont(42);
              tGlobal.SetTextAlign(13);
              tGlobal.SetTextSize(0.027);
              tGlobal.DrawLatex(0.045, 0.983, "Inclusive Fails (Preselection)");

              TPaveText* pCol1 = new TPaveText(0.34, 0.915, 0.52, 0.988, "NDC NB");
              pCol1->SetFillStyle(0);
              pCol1->SetBorderSize(0);
              pCol1->SetTextFont(42);
              pCol1->SetTextAlign(12);
              pCol1->SetTextSize(0.019);
              pCol1->AddText("A: #frac{E_{11}}{E_{33}} #geq 0.98");
              pCol1->AddText("B: et1 #leq 0.6");
              pCol1->Draw();
              keepTbl.push_back(pCol1);

              TPaveText* pCol2 = new TPaveText(0.51, 0.915, 0.67, 0.988, "NDC NB");
              pCol2->SetFillStyle(0);
              pCol2->SetBorderSize(0);
              pCol2->SetTextFont(42);
              pCol2->SetTextAlign(12);
              pCol2->SetTextSize(0.019);
              pCol2->AddText("C: et1 #geq 1.0");
              pCol2->AddText("D: et1 out of range");
              pCol2->Draw();
              keepTbl.push_back(pCol2);

              TPaveText* pCol3 = new TPaveText(0.66, 0.915, 0.82, 0.988, "NDC NB");
              pCol3->SetFillStyle(0);
              pCol3->SetBorderSize(0);
              pCol3->SetTextFont(42);
              pCol3->SetTextAlign(12);
              pCol3->SetTextSize(0.019);
              pCol3->AddText("E: #frac{E_{32}}{E_{35}} #leq 0.8");
              pCol3->AddText("F: #frac{E_{32}}{E_{35}} #geq 1.0");
              pCol3->Draw();
              keepTbl.push_back(pCol3);

              TPaveText* pCol4 = new TPaveText(0.80, 0.915, 0.975, 0.988, "NDC NB");
              pCol4->SetFillStyle(0);
              pCol4->SetBorderSize(0);
              pCol4->SetTextFont(42);
              pCol4->SetTextAlign(12);
              pCol4->SetTextSize(0.019);
              pCol4->AddText("G: #frac{E_{32}}{E_{35}} out of range");
              pCol4->AddText("H: w_{#eta}^{cogX} #geq 0.6");
              pCol4->Draw();
              keepTbl.push_back(pCol4);
              
              
              SaveCanvas(cTbl, JoinPath(outDir, "table3x3_preselectionFails.png"));

              for (auto* obj : keepTbl) delete obj;
          }

          // ---------------------------------------------------------------------------
          // additionalQA: survival-fraction vs pT (iso vs nonIso) + weta shape overlay
          // ---------------------------------------------------------------------------
          {
            const string qaDir = JoinPath(outDir, "additionalQA");
            EnsureDir(qaDir);

            // weta preselection cut (pass = weta < 0.6)
            const double kWetaCut  = 0.6;
            // et1 preselection window (pass = 0.6 < et1 < 1.0)
            const double kEt1Lo    = 0.6;
            const double kEt1Hi    = 1.0;

            const auto& bins = PtBins();

            // ------------------------------------------------------------------
            // Panel 1: survival fraction vs pT  (weta top, et1 bottom)
            // ------------------------------------------------------------------
            {
              vector<double> xPt(kNPtBins, 0.0), exPt(kNPtBins, 0.0);
              vector<double> fWIso(kNPtBins, 0.0),  fWNon(kNPtBins, 0.0);
              vector<double> fEIso(kNPtBins, 0.0),  fENon(kNPtBins, 0.0);

              for (int i = 0; i < kNPtBins; ++i)
              {
                const PtBin& b   = bins[i];
                xPt[i]  = 0.5 * (kPtEdges[(size_t)i] + kPtEdges[(size_t)i+1]);
                exPt[i] = 0.5 * (kPtEdges[(size_t)i+1] - kPtEdges[(size_t)i]);

                auto survFrac = [&](const string& hname, double lo, double hi, bool useLo)->double
                {
                  TH1* h = GetObj<TH1>(ds, hname, false, false, false);
                  if (!h) return -1.0;
                  const double tot = h->Integral(0, h->GetNbinsX()+1);
                  if (tot <= 0.0) return -1.0;
                  if (useLo)
                  {
                    // pass = val < lo
                    const int bCut = h->GetXaxis()->FindBin(lo - 1e-6);
                    return h->Integral(1, bCut) / tot;
                  }
                  else
                  {
                    // pass = lo < val < hi (open interval)
                    const int b1 = h->GetXaxis()->FindBin(lo + 1e-6);
                    const int b2 = h->GetXaxis()->FindBin(hi - 1e-6);
                    return h->Integral(b1, b2) / tot;
                  }
                };

                fWIso[i] = survFrac("h_ss_weta_iso"   + b.suffix, kWetaCut, 0.0,    true);
                fWNon[i] = survFrac("h_ss_weta_nonIso" + b.suffix, kWetaCut, 0.0,    true);
                fEIso[i] = survFrac("h_ss_et1_iso"    + b.suffix, kEt1Lo,   kEt1Hi, false);
                fENon[i] = survFrac("h_ss_et1_nonIso"  + b.suffix, kEt1Lo,   kEt1Hi, false);
              }

              // Build graphs (skip bins where hist was missing, flagged by -1)
              auto MakeGraph = [&](const vector<double>& yv,
                                   int col, int mStyle)->TGraphErrors*
              {
                vector<double> gx, gex, gy, gey;
                for (int i = 0; i < kNPtBins; ++i)
                {
                  if (yv[i] < 0.0) continue;
                  gx.push_back(xPt[i]);
                  gex.push_back(exPt[i]);
                  gy.push_back(yv[i]);
                  gey.push_back(0.0);
                }
                if (gx.empty()) return nullptr;
                TGraphErrors* g = new TGraphErrors((int)gx.size(),
                                                   &gx[0], &gy[0],
                                                   &gex[0], &gey[0]);
                g->SetLineColor(col);   g->SetLineWidth(2);
                g->SetMarkerColor(col); g->SetMarkerStyle(mStyle);
                g->SetMarkerSize(1.2);
                return g;
              };

              TGraphErrors* gWIso = MakeGraph(fWIso, kBlue+1,  20);
              TGraphErrors* gWNon = MakeGraph(fWNon, kRed+1,   21);
              TGraphErrors* gEIso = MakeGraph(fEIso, kBlue+1,  20);
              TGraphErrors* gENon = MakeGraph(fENon, kRed+1,   21);

              TCanvas cSurv("c_presel_survFrac","c_presel_survFrac", 1400, 600);
              cSurv.Divide(2,1, 0.002, 0.002);

                auto DrawSurvPad = [&](int padN, TGraphErrors* gIso, TGraphErrors* gNon,
                                       const char* varLabel, const char* cutText)
                {
                  cSurv.cd(padN);
                  gPad->SetLeftMargin(0.16);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.16);
                  gPad->SetTopMargin(0.10);
                  gPad->SetTicks(1,1);

                  // heap-allocate so the pad owns these past the lambda return
                  TH1F* hF = new TH1F(TString::Format("hSurvFrame_%d",padN),"",
                                      100, kPtEdges.front(), kPtEdges.back());
                  hF->SetDirectory(nullptr);
                  hF->SetStats(0);
                  hF->SetMinimum(0.0);
                  hF->SetMaximum(1.05);
                  hF->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  hF->GetYaxis()->SetTitle("Preselection survival fraction");
                  hF->GetXaxis()->SetTitleSize(0.055);
                  hF->GetYaxis()->SetTitleSize(0.052);
                  hF->GetXaxis()->SetLabelSize(0.048);
                  hF->GetYaxis()->SetLabelSize(0.048);
                  hF->Draw();

                  if (gIso) gIso->Draw("PE same");
                  if (gNon) gNon->Draw("PE same");

                  TLegend* leg = new TLegend(0.55, 0.20, 0.92, 0.38);
                  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.045);
                  if (gIso) leg->AddEntry(gIso, "Isolated",     "lpe");
                  if (gNon) leg->AddEntry(gNon, "Non-isolated", "lpe");
                  leg->Draw();

                  TLine lUnity(kPtEdges.front(), 1.0, kPtEdges.back(), 1.0);
                  lUnity.SetLineStyle(2); lUnity.SetLineColor(kGray+1); lUnity.SetLineWidth(1);
                  lUnity.DrawClone();

                  TLatex tt; tt.SetNDC(); tt.SetTextFont(42); tt.SetTextAlign(23); tt.SetTextSize(0.050);
                  tt.DrawLatex(0.50, 0.965, TString::Format("Presel survival: %s  (cut: %s)", varLabel, cutText));
                };

              DrawSurvPad(1, gWIso, gWNon, "weta", "weta < 0.6");
              DrawSurvPad(2, gEIso, gENon, "et1",  "0.6 < et1 < 1.0");

              SaveCanvas(cSurv, JoinPath(qaDir, "presel_survivalFraction_iso_vs_nonIso.png"));

              delete gWIso; delete gWNon; delete gEIso; delete gENon;
            }

            // ------------------------------------------------------------------
            // Panel 2: weta shape overlay (iso vs nonIso) for last two pT bins
            // ------------------------------------------------------------------
            {
              // last two bins: kNPtBins-2 and kNPtBins-1
              const int nPanel = std::min(2, kNPtBins);
              TCanvas cShape("c_weta_shape_lastBins","c_weta_shape_lastBins", 1400, 600);
              cShape.Divide(nPanel, 1, 0.002, 0.002);

              vector<TH1*> keepShape;
              keepShape.reserve(nPanel * 2);

              for (int ip = 0; ip < nPanel; ++ip)
              {
                const int i      = kNPtBins - nPanel + ip;
                const PtBin& b   = bins[i];

                TH1* hIso = GetObj<TH1>(ds, "h_ss_weta_iso"    + b.suffix, false, false, false);
                TH1* hNon = GetObj<TH1>(ds, "h_ss_weta_nonIso" + b.suffix, false, false, false);

                cShape.cd(ip + 1);
                gPad->SetLeftMargin(0.16);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.16);
                gPad->SetTopMargin(0.10);
                gPad->SetTicks(1,1);

                if (!hIso && !hNon)
                {
                  TLatex tm; tm.SetNDC(); tm.SetTextFont(42); tm.SetTextAlign(22); tm.SetTextSize(0.07);
                  tm.DrawLatex(0.50, 0.55, "MISSING");
                  continue;
                }

                TH1* cIso = hIso ? CloneTH1(hIso, TString::Format("wshp_iso_%d",i).Data())  : nullptr;
                TH1* cNon = hNon ? CloneTH1(hNon, TString::Format("wshp_non_%d",i).Data())  : nullptr;

                if (cIso) { cIso->SetDirectory(nullptr); cIso->Rebin(4); NormalizeToUnitArea(cIso);
                               cIso->SetLineColor(kBlue+1); cIso->SetLineWidth(2); keepShape.push_back(cIso); }
                if (cNon) { cNon->SetDirectory(nullptr); cNon->Rebin(4); NormalizeToUnitArea(cNon);
                               cNon->SetLineColor(kRed+1);  cNon->SetLineWidth(2); keepShape.push_back(cNon); }

                double ymax = 0.0;
                if (cIso) ymax = std::max(ymax, cIso->GetMaximum());
                if (cNon) ymax = std::max(ymax, cNon->GetMaximum());

                TH1* first = cIso ? cIso : cNon;
                first->SetTitle("");
                first->SetMaximum(ymax * 1.30);
                first->GetXaxis()->SetTitle("w_{#eta}^{cogX}");
                first->GetYaxis()->SetTitle("A.U.");
                first->GetXaxis()->SetTitleSize(0.055);
                first->GetYaxis()->SetTitleSize(0.052);
                first->GetXaxis()->SetLabelSize(0.048);
                first->GetYaxis()->SetLabelSize(0.048);
                first->Draw("hist");
                if (cIso && cNon) cNon->Draw("hist same");

                TLine lCut(kWetaCut, 0.0, kWetaCut, ymax * 1.25);
                lCut.SetLineColor(kBlack); lCut.SetLineStyle(2); lCut.SetLineWidth(2);
                lCut.DrawClone();

                TLegend leg(0.55, 0.68, 0.92, 0.84);
                leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.045);
                if (cIso) leg.AddEntry(cIso, "Isolated",     "l");
                if (cNon) leg.AddEntry(cNon, "Non-isolated", "l");
                leg.Draw();

                TLatex tt; tt.SetNDC(); tt.SetTextFont(42); tt.SetTextAlign(23); tt.SetTextSize(0.050);
                tt.DrawLatex(0.50, 0.965,
                  TString::Format("w_{#eta}^{cogX} shape  p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi));

                TLatex tc; tc.SetNDC(); tc.SetTextFont(42); tc.SetTextAlign(13); tc.SetTextSize(0.040);
                tc.DrawLatex(0.17, 0.86,
                  TString::Format("Dashed: presel cut (weta < %.1f)", kWetaCut));
              }

                SaveCanvas(cShape, JoinPath(qaDir, "weta_shape_iso_vs_nonIso_lastTwoBins.png"));
                for (auto* h : keepShape) delete h;
              }

              // ------------------------------------------------------------------
              // Panel 3: ABCD discrimination diagnostics vs pT
              //   (a) f_tight_iso = A/(A+C)  vs  f_tight_nonIso = B/(B+D)
              //       If these converge => tight-ID loses discrimination => purity->0
              //   (b) R = A*D/(B*C): ABCD signal ratio; R=1 => purity=0
              // All computed purely from existing ABCD count histograms.
              // ------------------------------------------------------------------
              {
                const double kWetaCutABCD = 0.6;

                  const char* regionNames[4]    = {"A (iso&tight)",    "B (nonIso&tight)",
                                                   "C (iso&nonTight)", "D (nonIso&nonTight)"};
                    vector<double> xPtD(kNPtBins, 0.0), exPtD(kNPtBins, 0.0);
                    vector<double> yFtIso(kNPtBins, -1.0), eyFtIso(kNPtBins, 0.0);
                    vector<double> yFtNon(kNPtBins, -1.0), eyFtNon(kNPtBins, 0.0);
                    vector<double> yR(kNPtBins,     -1.0), eyR(kNPtBins,     0.0);

                    for (int i = 0; i < kNPtBins; ++i)
                    {
                      const PtBin& b = bins[i];
                      xPtD[i]  = 0.5 * (kPtEdges[(size_t)i] + kPtEdges[(size_t)i+1]);
                      exPtD[i] = 0.5 * (kPtEdges[(size_t)i+1] - kPtEdges[(size_t)i]);

                      const double A  = Read1BinCount(ds, "h_isIsolated_isTight"   + b.suffix);
                      const double B  = Read1BinCount(ds, "h_notIsolated_isTight"  + b.suffix);
                      const double C  = Read1BinCount(ds, "h_isIsolated_notTight"  + b.suffix);
                      const double D  = Read1BinCount(ds, "h_notIsolated_notTight" + b.suffix);

                      // f_tight,iso = A/(A+C)
                      // binomial error: sigma = sqrt(A*C) / (A+C)^(3/2)
                      if ((A + C) > 0.0)
                      {
                        yFtIso[i]  = A / (A + C);
                        eyFtIso[i] = (A > 0.0 && C > 0.0)
                                     ? std::sqrt(A * C) / std::pow(A + C, 1.5)
                                     : std::sqrt(A) / (A + C);
                      }

                      // f_tight,nonIso = B/(B+D)
                      if ((B + D) > 0.0)
                      {
                        yFtNon[i]  = B / (B + D);
                        eyFtNon[i] = (B > 0.0 && D > 0.0)
                                     ? std::sqrt(B * D) / std::pow(B + D, 1.5)
                                     : std::sqrt(B) / (B + D);
                      }

                      // R = A*D/(B*C)
                      // relative error: sigma_R/R = sqrt(1/A + 1/B + 1/C + 1/D)
                      if (A > 0.0 && B > 0.0 && C > 0.0 && D > 0.0)
                      {
                        yR[i]  = (A * D) / (B * C);
                        eyR[i] = yR[i] * std::sqrt(1.0/A + 1.0/B + 1.0/C + 1.0/D);
                      }
                    }

                    auto MakeGrE = [&](const vector<double>& yv, const vector<double>& eyv,
                                       int col, int ms)->TGraphErrors*
                    {
                      vector<double> gx, gex, gy, gey;
                      for (int i = 0; i < kNPtBins; ++i)
                      {
                        if (yv[i] < 0.0) continue;
                        gx.push_back(xPtD[i]); gex.push_back(exPtD[i]);
                        gy.push_back(yv[i]);   gey.push_back(eyv[i]);
                      }
                      if (gx.empty()) return nullptr;
                      TGraphErrors* g = new TGraphErrors((int)gx.size(),
                                                         &gx[0], &gy[0], &gex[0], &gey[0]);
                      g->SetLineColor(col); g->SetLineWidth(2);
                      g->SetMarkerColor(col); g->SetMarkerStyle(ms); g->SetMarkerSize(1.3);
                      return g;
                    };

                    TGraphErrors* gFtIso = MakeGrE(yFtIso, eyFtIso, kBlue+1, 20);
                    TGraphErrors* gFtNon = MakeGrE(yFtNon, eyFtNon, kRed+1,  21);
                    TGraphErrors* gR     = MakeGrE(yR,     eyR,     kBlack,  20);

                    // --- canvas 1: tight fraction iso vs nonIso ---
                    {
                      TCanvas cDisc1("c_abcd_tightFrac","c_abcd_tightFrac", 900, 700);
                      cDisc1.SetLeftMargin(0.16); cDisc1.SetRightMargin(0.05);
                      cDisc1.SetBottomMargin(0.14); cDisc1.SetTopMargin(0.10); cDisc1.SetTicks(1,1);

                      TH1F* hF1 = new TH1F("hDiscFrame1","", 100, kPtEdges.front(), kPtEdges.back());
                      hF1->SetDirectory(nullptr); hF1->SetStats(0);
                      hF1->SetMinimum(0.0); hF1->SetMaximum(1.05);
                      hF1->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hF1->GetYaxis()->SetTitle("Tight fraction");
                      hF1->GetXaxis()->SetTitleSize(0.052); hF1->GetYaxis()->SetTitleSize(0.050);
                      hF1->GetXaxis()->SetLabelSize(0.046); hF1->GetYaxis()->SetLabelSize(0.046);
                      hF1->Draw();

                        TLine* lHalfPtr = new TLine(kPtEdges.front(), 0.5, kPtEdges.back(), 0.5);
                        lHalfPtr->SetLineStyle(2); lHalfPtr->SetLineColor(kGray+1); lHalfPtr->SetLineWidth(2);
                        lHalfPtr->Draw();

                        if (gFtIso) gFtIso->Draw("PE same");
                        if (gFtNon) gFtNon->Draw("PE same");

                        {
                          TLegend* leg1 = new TLegend(0.18, 0.74, 0.62, 0.89);
                          leg1->SetBorderSize(0); leg1->SetFillStyle(0);
                          leg1->SetTextFont(42);  leg1->SetTextSize(0.036);
                          if (gFtIso) leg1->AddEntry(gFtIso, "f_{tight,iso}  = A/(A+C)", "lpe");
                          if (gFtNon) leg1->AddEntry(gFtNon, "f_{tight,nonIso} = B/(B+D)", "lpe");
                          leg1->AddEntry(lHalfPtr, "No-discrimination boundary (f = 0.5)", "l");
                          leg1->Draw();
                        }

                      {
                        TLatex tt1; tt1.SetNDC(); tt1.SetTextFont(42); tt1.SetTextAlign(23); tt1.SetTextSize(0.048);
                        tt1.DrawLatex(0.50, 0.965, "Tight-ID discrimination vs p_{T}^{#gamma}");
                      }

                      SaveCanvas(cDisc1, JoinPath(qaDir, "abcd_tightFrac_discrimination_vs_pT.png"));
                    }

                    // --- canvas 2: ABCD signal ratio R = A*D/(B*C) ---
                    {
                      double Rmax = 1.5;
                      for (int i = 0; i < kNPtBins; ++i)
                        if (yR[i] > 0.0) Rmax = std::max(Rmax, (yR[i] + eyR[i]) * 1.25);
                      Rmax = std::min(Rmax, 15.0);

                      TCanvas cDisc2("c_abcd_Rratio","c_abcd_Rratio", 900, 700);
                      cDisc2.SetLeftMargin(0.16); cDisc2.SetRightMargin(0.05);
                      cDisc2.SetBottomMargin(0.14); cDisc2.SetTopMargin(0.10); cDisc2.SetTicks(1,1);

                      TH1F* hF2 = new TH1F("hDiscFrame2","", 100, kPtEdges.front(), kPtEdges.back());
                      hF2->SetDirectory(nullptr); hF2->SetStats(0);
                      hF2->SetMinimum(0.0); hF2->SetMaximum(Rmax);
                      hF2->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hF2->GetYaxis()->SetTitle("R = A#timesD / (B#timesC)");
                      hF2->GetXaxis()->SetTitleSize(0.052); hF2->GetYaxis()->SetTitleSize(0.050);
                      hF2->GetXaxis()->SetLabelSize(0.046); hF2->GetYaxis()->SetLabelSize(0.046);
                      hF2->Draw();

                      {
                        TLine lOne(kPtEdges.front(), 1.0, kPtEdges.back(), 1.0);
                        lOne.SetLineStyle(2); lOne.SetLineColor(kRed+1); lOne.SetLineWidth(2);
                        lOne.DrawClone();
                      }
                      if (gR) gR->Draw("PE same");

                      {
                        TLatex tt2; tt2.SetNDC(); tt2.SetTextFont(42); tt2.SetTextAlign(23); tt2.SetTextSize(0.048);
                        tt2.DrawLatex(0.50, 0.965, "ABCD signal ratio R = A#timesD/(B#timesC)");
                        TLatex tc2; tc2.SetNDC(); tc2.SetTextFont(42); tc2.SetTextAlign(13); tc2.SetTextSize(0.038);
                        tc2.DrawLatex(0.17, 0.88, "Purity = 1 #minus 1/R  (R=1 #Rightarrow purity=0)");
                        tc2.DrawLatex(0.17, 0.82, "Dashed red: R=1 boundary");
                      }

                      SaveCanvas(cDisc2, JoinPath(qaDir, "abcd_signalRatio_R_vs_pT.png"));
                    }

                    delete gFtIso; delete gFtNon; delete gR;
              }
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

        const bool doEmbeddedVariantOverlay = (isSimEmbedded && ds.isSim && !ds.centFolder.empty());

        string outDir;
        if (doEmbeddedVariantOverlay) outDir = JoinPath(ds.outBase, "isoQA");
        else if (ds.isSim)           outDir = JoinPath(ds.outBase, "isoQAgeneral");
        else                         outDir = JoinPath(ds.outBase, "baselineData/isoQAgeneral");

        EnsureDir(outDir);
        for (const auto& b : PtBins()) EnsureDir(JoinPath(outDir, b.folder));

        // ------------------------------------------------------------------
        // Embedded photon20 special mode:
        //   ONLY produce reco total-isolation overlays across UE variants
        //   inside embeddedPhoton20/<cent>/isoQA/
        // ------------------------------------------------------------------
        if (doEmbeddedVariantOverlay)
        {
            struct EmbeddedIsoVariantHandle
            {
              string variant;
              int color = 1;
              TFile* file = nullptr;
              Dataset dsVar;
            };

            vector<EmbeddedIsoVariantHandle> vars;
            vars.reserve(4);

            auto addVariant = [&](const string& variant, int color)
            {
              EmbeddedIsoVariantHandle H;
              H.variant = variant;
              H.color = color;
              H.file = TFile::Open(InputSimEmbedded(variant).c_str(), "READ");

                if (!H.file || H.file->IsZombie())
                {
                  if (H.file) { H.file->Close(); delete H.file; H.file = nullptr; }

                  cout << ANSI_BOLD_YEL
                       << "[WARN] Missing embedded iso overlay input: "
                       << InputSimEmbedded(variant)
                       << ANSI_RESET << "\n";

                  vars.push_back(std::move(H));
                  return;
                }

                H.dsVar.label = "SIM";
                H.dsVar.isSim = true;
                H.dsVar.trigger = "";
                H.dsVar.topDirName = kDirSIM;
                H.dsVar.inFilePath = InputSimEmbedded(variant);
                H.dsVar.outBase = outDir;
                H.dsVar.centFolder = ds.centFolder;
                H.dsVar.centSuffix = ds.centSuffix;
                H.dsVar.centLabel = ds.centLabel;
                H.dsVar.file = H.file;
                H.dsVar.topDir = H.file->GetDirectory(H.dsVar.topDirName.c_str());

                if (!H.dsVar.topDir)
                {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] Missing " << kDirSIM << " directory in embedded iso overlay input: "
                       << InputSimEmbedded(variant)
                     << ANSI_RESET << "\n";
              }

              vars.push_back(std::move(H));
            };

            addVariant("noSub",       kBlack);
            addVariant("baseVariant", kBlue+1);
            addVariant("variantA",    kRed+1);
            addVariant("variantB",    kGreen+2);

            auto DrawVariantOverlayPad =
              [&](const PtBin& b, bool compact, vector<TObject*>* keepAlive)
            {
              vector<TH1*> hs;
              vector<string> labels;
              double ymax = 0.0;

              for (auto& V : vars)
              {
                if (!V.dsVar.topDir) continue;

                TH1* hSrc = GetObj<TH1>(V.dsVar, "h_Eiso" + b.suffix, false, false, false);
                if (!hSrc) continue;

                TH1* h = CloneTH1(
                  hSrc,
                  TString::Format("h_Eiso_variantOverlay_%s_%s_%s",
                    ds.centFolder.c_str(), b.folder.c_str(), V.variant.c_str()).Data()
                );
                if (!h) continue;

                h->SetLineColor(V.color);
                h->SetLineWidth(2);
                h->SetFillStyle(0);

                ymax = std::max(ymax, h->GetMaximum());
                hs.push_back(h);
                labels.push_back(V.variant);

                if (keepAlive) keepAlive->push_back(h);
              }

              if (hs.empty())
              {
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextSize(compact ? 0.060 : 0.055);
                t.DrawLatex(0.18, 0.55, "MISSING");
                return;
              }

              TH1* first = hs[0];
              first->SetTitle("");
              first->GetXaxis()->SetTitle("E_{iso} [GeV]");
              first->GetYaxis()->SetTitle("Counts");
              first->GetXaxis()->SetTitleSize(compact ? 0.060 : 0.055);
              first->GetYaxis()->SetTitleSize(compact ? 0.060 : 0.055);
              first->GetXaxis()->SetLabelSize(compact ? 0.045 : 0.045);
              first->GetYaxis()->SetLabelSize(compact ? 0.045 : 0.045);
              first->GetYaxis()->SetTitleOffset(compact ? 1.00 : 1.15);
              first->SetMinimum(0.0);
              first->SetMaximum((ymax > 0.0) ? (1.25 * ymax) : 1.0);

              first->Draw("hist");
              for (std::size_t ih = 1; ih < hs.size(); ++ih) hs[ih]->Draw("hist same");

              TLegend* leg = new TLegend(compact ? 0.52 : 0.56, compact ? 0.63 : 0.68, 0.92, 0.90);
              leg->SetBorderSize(0);
              leg->SetFillStyle(0);
              leg->SetTextFont(42);
              leg->SetTextSize(compact ? 0.038 : 0.032);
              for (std::size_t ih = 0; ih < hs.size(); ++ih)
              {
                leg->AddEntry(hs[ih], labels[ih].c_str(), "l");
              }
              leg->Draw();
              if (keepAlive) keepAlive->push_back(leg);

              if (compact)
              {
                TLatex tt;
                tt.SetTextFont(42);
                tt.SetNDC();
                tt.SetTextAlign(23);
                tt.SetTextSize(0.052);
                tt.DrawLatex(0.50, 0.965,
                  TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
              }
              else
              {
                vector<string> lines = DefaultHeaderLines(ds);
                lines.push_back("E_{T}^{iso, Total} overlay by embedded UE variant");
                lines.push_back(TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                DrawLatexLines(0.16, 0.92, lines, 0.030, 0.040);
              }
            };

            for (int i = 0; i < kNPtBins; ++i)
            {
              const PtBin& b = PtBins()[i];

              TCanvas c(
                TString::Format("c_embIsoVariantOverlay_%s_%s", ds.centFolder.c_str(), b.folder.c_str()).Data(),
                "c_embIsoVariantOverlay", 900, 700
              );
              ApplyCanvasMargins1D(c);
              c.cd();

              vector<TObject*> keep;
              DrawVariantOverlayPad(b, false, &keep);

              const string fp = JoinPath(outDir, b.folder + "/Eiso_total_variantOverlay_" + b.folder + ".png");
              SaveCanvas(c, fp);

              for (auto* obj : keep) delete obj;
            }

            const int nCols   = 3;
            const int nRows   = 3;
            const int perPage = nCols * nRows;
            const int nPages  = (kNPtBins + perPage - 1) / perPage;

            for (int ipage = 0; ipage < nPages; ++ipage)
            {
              TCanvas cTbl(
                TString::Format("c_embIsoVariantOverlayTbl_%s_page%d", ds.centFolder.c_str(), ipage + 1).Data(),
                "c_embIsoVariantOverlayTbl", 1800, 1400
              );
              cTbl.Divide(nCols, nRows, 0.001, 0.001);

              vector<TObject*> keepTbl;

              for (int ipad = 0; ipad < perPage; ++ipad)
              {
                cTbl.cd(ipad + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.10);
                gPad->SetTicks(1,1);

                const int idx = ipage * perPage + ipad;
                if (idx >= kNPtBins)
                {
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.060);
                  t.DrawLatex(0.32, 0.55, "EMPTY");
                  continue;
                }

                DrawVariantOverlayPad(PtBins()[idx], true, &keepTbl);
              }

              SaveCanvas(
                cTbl,
                JoinPath(outDir,
                  TString::Format("table3x3_Eiso_total_variantOverlay_page%d.png", ipage + 1).Data())
              );

              for (auto* obj : keepTbl) delete obj;
            }

            for (auto& V : vars)
            {
              if (V.file)
              {
                V.file->Close();
                delete V.file;
                V.file = nullptr;
                V.dsVar.file = nullptr;
                V.dsVar.topDir = nullptr;
              }
            }

            return;
        }

        // 1) RECO/DATA isolation QA (always runs outside embedded special mode)
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

            Make3x3Table_TH1(ds,
                             "h_nIsoTightPhoCand",
                             outDirMult,
                             "table3x3_nIsoTightPhoCand.png",
                             "N_{#gamma}^{iso+tight} candidates",
                             "A.U.",
                             false,
                             true,
                             multLines);

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

                const int b2 = h->GetXaxis()->FindBin(2.0);
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

        TCanvas c(
          TString::Format("c_abcd_cnt_tbl_%s", ds.label.c_str()).Data(),
          "c_abcd_cnt_tbl", 2700, 2300
        );
        c.Divide(3,3, 0.001, 0.001);

        vector<TObject*> keepAlive;
        keepAlive.reserve(kNPtBins * (1 + 4));

        const auto& bins = PtBins();

        for (int i = 0; i < kNPtBins; ++i)
        {
          c.cd(i+1);
          if (!gPad) continue;

          gPad->SetLeftMargin(0.16);
          gPad->SetRightMargin(0.05);
          gPad->SetTopMargin(0.11);
          gPad->SetBottomMargin(0.22);
          gPad->SetTicks(1,1);

          if (i < 3)
          {
            const double x1 = gPad->GetXlowNDC();
            const double y1 = gPad->GetYlowNDC();
            const double x2 = x1 + gPad->GetWNDC();
            const double y2 = y1 + gPad->GetHNDC();
            const double shiftDown = 0.045;
            gPad->SetPad(x1, y1 - shiftDown, x2, y2 - shiftDown);
          }

          const PtBin& b   = bins[i];
          const string suf = b.suffix;

          const double A = Read1BinCount(ds, "h_isIsolated_isTight"     + suf);
          const double B = Read1BinCount(ds, "h_notIsolated_isTight"    + suf);
          const double Cc = Read1BinCount(ds, "h_isIsolated_notTight"   + suf);
          const double D = Read1BinCount(ds, "h_notIsolated_notTight"   + suf);

          const double vals[4] = {A, B, Cc, D};

          double ymax = 0.0;
          for (int ib = 0; ib < 4; ++ib) ymax = std::max(ymax, vals[ib]);
          const double yMaxPlot = (ymax > 0.0) ? (1.35 * ymax) : 1.0;

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
          hAxis->GetXaxis()->LabelsOption("h");
          hAxis->GetXaxis()->SetLabelSize(0.070);
          hAxis->GetXaxis()->SetLabelOffset(0.010);

          hAxis->GetYaxis()->SetTitleSize(0.054);
          hAxis->GetYaxis()->SetTitleOffset(1.28);
          hAxis->GetYaxis()->SetLabelSize(0.046);

          hAxis->SetLineColor(1);
          hAxis->SetLineWidth(2);
          hAxis->SetFillStyle(0);
          hAxis->Draw("hist");

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

          TLatex t;
          t.SetTextFont(42);
          t.SetTextAlign(22);
          t.SetTextSize(0.052);

          for (int ib = 1; ib <= 4; ++ib)
          {
            const double y = vals[ib-1];
            if (y <= 0.0) continue;

            const double x = hAxis->GetXaxis()->GetBinCenter(ib);
            const double yText = std::min(y + 0.025*yMaxPlot, 0.93*yMaxPlot);
            t.DrawLatex(x, yText, TString::Format("%.0f", y).Data());
          }

          TLatex tt;
          tt.SetTextFont(42);
          tt.SetNDC();
          tt.SetTextAlign(23);
          tt.SetTextSize(0.060);
          tt.DrawLatex(0.50, 0.965,
              TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data()
          );

          keepAlive.push_back(hAxis);
        }

        c.cd(0);

        TLatex tGlobal;
        tGlobal.SetNDC();
        tGlobal.SetTextFont(42);
        tGlobal.SetTextAlign(22);
        tGlobal.SetTextSize(0.028);
        tGlobal.DrawLatex(0.50, 0.965,
          "ABCD counts after preselection (A=iso&tight, B=nonIso&tight, C=iso&nonTight, D=nonIso&nonTight)");

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

          // purity plots: use physical p_{T}^{#gamma} on the x axis and explicit purity uncertainties
          {
            auto RawPurityError = [&](double A, double B, double C, double D)->double
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
            };

            auto CorrPurityValue = [&](double A, double B, double C, double D, int idx)->double
            {
              double Praw = 0.0;
              if (A > 0.0 && D > 0.0)
              {
                double Asig = A - B * (C / D);
                if (Asig < 0.0) Asig = 0.0;
                Praw = Asig / A;
              }

              if (!lf.available) return Praw;

              double SA = 0.0;
              const bool ok = SolveLeakageCorrectedSA(A, B, C, D, lf.fB[idx], lf.fC[idx], lf.fD[idx], SA);
              if (ok && A > 0.0) return SA / A;

              return Praw;
            };

            auto CorrPurityError = [&](double A, double B, double C, double D, int idx)->double
            {
              if (!(lf.available && hasCorr[idx])) return RawPurityError(A, B, C, D);
              if (A <= 0.0) return 0.0;

              const double dA = std::sqrt(std::max(A, 1.0));
              const double dB = std::sqrt(std::max(B, 1.0));
              const double dC = std::sqrt(std::max(C, 1.0));
              const double dD = std::sqrt(std::max(D, 1.0));

              const double Aup = A + dA;
              const double Adn = std::max(0.0, A - dA);
              const double Bup = B + dB;
              const double Bdn = std::max(0.0, B - dB);
              const double Cup = C + dC;
              const double Cdn = std::max(0.0, C - dC);
              const double Dup = D + dD;
              const double Ddn = std::max(0.0, D - dD);

              const double dPdA = (Aup > Adn)
                ? (CorrPurityValue(Aup, B,   C,   D,   idx) - CorrPurityValue(Adn, B,   C,   D,   idx)) / (Aup - Adn)
                : 0.0;
              const double dPdB = (Bup > Bdn)
                ? (CorrPurityValue(A,   Bup, C,   D,   idx) - CorrPurityValue(A,   Bdn, C,   D,   idx)) / (Bup - Bdn)
                : 0.0;
              const double dPdC = (Cup > Cdn)
                ? (CorrPurityValue(A,   B,   Cup, D,   idx) - CorrPurityValue(A,   B,   Cdn, D,   idx)) / (Cup - Cdn)
                : 0.0;
              const double dPdD = (Dup > Ddn)
                ? (CorrPurityValue(A,   B,   C,   Dup, idx) - CorrPurityValue(A,   B,   C,   Ddn, idx)) / (Dup - Ddn)
                : 0.0;

              double var = 0.0;
              if (A > 0.0) var += dPdA * dPdA * A;
              if (B > 0.0) var += dPdB * dPdB * B;
              if (C > 0.0) var += dPdC * dPdC * C;
              if (D > 0.0) var += dPdD * dPdD * D;

              return (var > 0.0) ? std::sqrt(var) : 0.0;
            };

            vector<double> x(kNPtBins, 0.0), ex(kNPtBins, 0.0);
            vector<double> yRaw(kNPtBins, 0.0), eyRaw(kNPtBins, 0.0);
            vector<double> yCorr(kNPtBins, 0.0), eyCorr(kNPtBins, 0.0);

            for (int i = 0; i < kNPtBins; ++i)
            {
              const PtBin& b = PtBins()[i];
              const string& suf = b.suffix;

              const double A = Read1BinCount(ds, "h_isIsolated_isTight" + suf);
              const double B = Read1BinCount(ds, "h_notIsolated_isTight" + suf);
              const double C = Read1BinCount(ds, "h_isIsolated_notTight" + suf);
              const double D = Read1BinCount(ds, "h_notIsolated_notTight" + suf);

              const double ptLo = kPtEdges[(std::size_t)i];
              const double ptHi = kPtEdges[(std::size_t)i + 1];

              x[i]  = 0.5 * (ptLo + ptHi);
              ex[i] = 0.5 * (ptHi - ptLo);

              yRaw[i]  = purityRaw[i];
              eyRaw[i] = RawPurityError(A, B, C, D);

              yCorr[i]  = purityCorr[i];
              eyCorr[i] = CorrPurityError(A, B, C, D, i);
            }

            // purity_raw plot
            {
              TGraphErrors gRaw(kNPtBins, &x[0], &yRaw[0], &ex[0], &eyRaw[0]);
              gRaw.SetLineWidth(2);
              gRaw.SetLineColor(kBlack);
              gRaw.SetMarkerStyle(20);
              gRaw.SetMarkerSize(1.2);
              gRaw.SetMarkerColor(kBlack);

              TCanvas c("c_pur_raw","c_pur_raw",900,700);
              ApplyCanvasMargins1D(c);

              TH1F hFrame("hPurRawFrame","",100, kPtEdges.front(), kPtEdges.back());
              hFrame.SetDirectory(nullptr);
              hFrame.SetStats(0);
              hFrame.SetMinimum(0.0);
              hFrame.SetMaximum(1.05);
              hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hFrame.GetYaxis()->SetTitle("Purity (raw ABCD)");
              hFrame.Draw();

              gRaw.Draw("P SAME");

              vector<string> box;
              box.push_back("ABCD purity (raw)");
              box.push_back("Cuts: p_{T}^{#gamma} #geq 5 GeV, |#eta| < 0.7, preselection pass");
              box.push_back("Iso: E_{iso} < 1.08128 + 0.0299107 p_{T}^{#gamma}");
              box.push_back("NonIso: E_{iso} > isoThresh + 1 GeV");
              DrawLatexLines(0.2, 0.92, DefaultHeaderLines(ds), 0.034, 0.045);
              DrawLatexLines(0.2, 0.8, box, 0.030, 0.040);

              const string fp = JoinPath(outDir, ds.isSim ? "purity_raw_SIM.png" : "purity_raw_DATA.png");
              SaveCanvas(c, fp);
            }

            // overlay raw vs corrected if correction used
            bool anyCorr = false;
            for (bool b : hasCorr) if (b) { anyCorr = true; break; }

            if (anyCorr)
            {
              TGraphErrors gRaw(kNPtBins, &x[0], &yRaw[0],  &ex[0], &eyRaw[0]);
              TGraphErrors gCor(kNPtBins, &x[0], &yCorr[0], &ex[0], &eyCorr[0]);

              gRaw.SetLineWidth(2);
              gRaw.SetLineColor(kBlack);
              gRaw.SetMarkerStyle(20);
              gRaw.SetMarkerSize(1.2);
              gRaw.SetMarkerColor(kBlack);

              gCor.SetLineWidth(2);
              gCor.SetLineColor(kBlue + 1);
              gCor.SetMarkerStyle(24);
              gCor.SetMarkerSize(1.2);
              gCor.SetMarkerColor(kBlue + 1);

              TCanvas c("c_pur_ov","c_pur_ov",900,700);
              ApplyCanvasMargins1D(c);

              TH1F hFrame("hPurCorFrame","",100, kPtEdges.front(), kPtEdges.back());
              hFrame.SetDirectory(nullptr);
              hFrame.SetStats(0);
              hFrame.SetMinimum(0.0);
              hFrame.SetMaximum(1.05);
              hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hFrame.GetYaxis()->SetTitle("Purity");
              hFrame.Draw();

              gRaw.Draw("P SAME");
              gCor.Draw("P SAME");

              TLegend leg(0.62,0.77,0.92,0.90);
              leg.SetTextFont(42);
              leg.SetTextSize(0.033);
              leg.AddEntry(&gRaw, "Raw ABCD", "lp");
              leg.AddEntry(&gCor, "Leakage-corrected", "lp");
              leg.Draw();

              const double etaCut = kPhotonEtaAbsMax;

              TLatex tTitle;
              tTitle.SetNDC(true);
              tTitle.SetTextFont(42);
              tTitle.SetTextAlign(23);
              tTitle.SetTextSize(0.045);
              tTitle.DrawLatex(0.50, 0.96,
                  ds.isSim
                    ? "Purity for SIM"
                    : "Purity for Run24pp, Photon 4 GeV + MBD NS #geq 1");

              TLatex tCuts;
              tCuts.SetNDC(true);
              tCuts.SetTextFont(42);
              tCuts.SetTextAlign(13);
              tCuts.SetTextSize(0.038);
              tCuts.DrawLatex(0.162, 0.88, TString::Format("|v_{z}| < %.3g cm", std::fabs(vzCutCm)).Data());
              tCuts.DrawLatex(0.162, 0.82, TString::Format("|#eta^{#gamma}| < %.3g", etaCut).Data());

              const string fp = JoinPath(outDir, ds.isSim
                  ? "purity_raw_vs_leakageCorrected_SIM.png"
                  : "purity_raw_vs_leakageCorrected_DATA.png");

              SaveCanvas(c, fp);
            }
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

                      std::vector<std::string> cutLines;
                      cutLines.push_back(TString::Format("Back-to-back: #Delta#phi_{#gamma,jet} > %s", B2BLabel().c_str()).Data());
                      cutLines.push_back(TString::Format("p_{T}^{jet} > %.0f GeV", static_cast<double>(kJetPtMin)).Data());
                      DrawLatexLines(0.70, 0.74, cutLines, 0.030, 0.038);

                      SaveCanvas(can, JoinPath(outDirHere, TString::Format("overlay_pTbin%d.png", ib).Data()));

                      if (a) delete a;
                      if (b) delete b;
                      if (c6) delete c6;
                  }

                  // Single 3x3 summary table using the full 9 JES3 pT bins
                  const int nCols   = 3;
                  const int nRows   = 3;
                  const int perPage = nCols * nRows; // 9

                  const int startBinForTable = 1;
                  const int nTableBins = std::min(perPage, nPt);

                  if (nTableBins > 0)
                  {
                      TCanvas canTbl(
                        TString::Format("c_tbl_%s_all9", outDirHere.c_str()).Data(),
                        "c_tbl_overlay_all9", 1500, 1200
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

                          TLatex tcuts;
                          tcuts.SetNDC(true);
                          tcuts.SetTextFont(42);
                          tcuts.SetTextAlign(33);
                          tcuts.SetTextSize(0.040);
                          tcuts.DrawLatex(0.92, 0.62, TString::Format("|#Delta#phi(#gamma,jet)| > %s", B2BLabel().c_str()).Data());
                          tcuts.DrawLatex(0.92, 0.54, TString::Format("Reco (p_{T}^{jet} > %.0f GeV)", static_cast<double>(kJetPtMin)).Data());

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
                        ? TString::Format("table3x3_overlay_alphaLT%.2f.png", alphaMax).Data()
                        : "table3x3_overlay_integratedAlpha.png";

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

                  EnsureMerged(MergedSimPath("photonJet5and10merged_SIM", "RecoilJets_photonjet5plus10_MERGED.root"),
                                                 {InputSim("photonjet5"), InputSim("photonjet10")},
                                                 {kSigmaPhoton5_pb, kSigmaPhoton10_pb},
                                                 {"photonJet5", "photonJet10"});

                  EnsureMerged(MergedSimPath("photonJet5and20merged_SIM", "RecoilJets_photonjet5plus20_MERGED.root"),
                             {InputSim("photonjet5"), InputSim("photonjet20")},
                             {kSigmaPhoton5_pb, kSigmaPhoton20_pb},
                             {"photonJet5", "photonJet20"});

                  EnsureMerged(MergedSimPath("photonJet10and20merged_SIM", "RecoilJets_photonjet10plus20_MERGED.root"),
                             {InputSim("photonjet10"), InputSim("photonjet20")},
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
                      {"photon jet 5 GeV",            InputSim("photonjet5"),                                                                          kBlack,      false, nullptr, nullptr},
                      {"photon jet 10 GeV",           InputSim("photonjet10"),                                                                         kRed+1,      false, nullptr, nullptr},
                      {"photon jet 20 GeV",           InputSim("photonjet20"),                                                                         kBlue+1,     false, nullptr, nullptr},
                      {"photon jet 5 + 10 GeV",       MergedSimPath("photonJet5and10merged_SIM", "RecoilJets_photonjet5plus10_MERGED.root"),            kGreen+3,    true,  nullptr, nullptr},
                      {"photon jet 5 + 20 GeV",       MergedSimPath("photonJet5and20merged_SIM", "RecoilJets_photonjet5plus20_MERGED.root"),            kMagenta+1,  true,  nullptr, nullptr},
                      {"photon jet 10 + 20 GeV",      MergedSimPath("photonJet10and20merged_SIM", "RecoilJets_photonjet10plus20_MERGED.root"),          kOrange+7,   true,  nullptr, nullptr},
                      {"photon jet 5 + 10 + 20 GeV",  MergedSimPath("photonJet5and10and20merged_SIM", "RecoilJets_photonjet5plus10plus20_MERGED.root"), kViolet+1,   true,  nullptr, nullptr}
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
              JoinPath(outRecoDir, "table3x3_overlay_integratedAlpha_dPhiCuts.png");

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
            const std::string altKey = CfgTag();
            const std::string altIn10 = InputSim("photonjet10");
            const std::string altIn20 = InputSim("photonjet20");
            const std::string altMerged = MergedSimPath("photonJet10and20merged_SIM", "RecoilJets_photonjet10plus20_MERGED.root");

            const double jetMinPtGeV = static_cast<double>(kJetPtMin);
            const std::string bbLabel = B2BLabel();

            auto EnsureAltMerged = [&]()->bool
            {
              if (!gSystem->AccessPathName(altMerged.c_str())) return true; // exists

              cout << ANSI_BOLD_CYN
                   << "\n[MERGE] Building merged SIM10+20 file for Δφ overlay:\n"
                   << "  cfgKey   = " << altKey << "\n"
                   << "  in10     = " << altIn10 << "\n"
                   << "  in20     = " << altIn20 << "\n"
                   << "  out      = " << altMerged << "\n"
                   << ANSI_RESET;

              return BuildMergedSIMFile_PhotonSlices(
                {altIn10, altIn20},
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
            // Make a 3x3 table (9 pads) across pTgamma bins: ProjectY (xJ) integrating alpha
            // AND also write per-pT individual overlays (same styling) into perPtBin_dPhiCuts/
            // -------------------------------------------------------------------------
            const int nCols = 3;
            const int nRows = 3;

            const double R = RFromKey(rKey);

            const int nBase = hBase->GetXaxis()->GetNbins();
            const int nAlt  = hAlt->GetXaxis()->GetNbins();
            const int nPt   = std::min(nBase, nAlt);

            const std::string outPerPtDir = JoinPath(outRecoDir, "perPtBin_dPhiCuts");
            EnsureDir(outPerPtDir);

            TCanvas c("c_tbl_dphiCuts_r04", "c_tbl_dphiCuts_r04", 1500, 1200);
            c.Divide(nCols, nRows, 0.001, 0.001);

            std::vector<TH1*> keep;
            keep.reserve(2 * 9);

            const int ibStart = 1;
            const int nAvail  = (nPt >= ibStart) ? (nPt - ibStart + 1) : 0;
            const int nPads   = std::min(9, nAvail);

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
                if (hPi2) leg.AddEntry(hPi2,  TString::Format("|#Delta#phi(#gamma,jet)| > %s", baseLabel.c_str()).Data(), "ep");
                if (h7p8) leg.AddEntry(h7p8,  TString::Format("|#Delta#phi(#gamma,jet)| > %s", bbLabel.c_str()).Data(), "ep");
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
                  JoinPath(outRecoDir, "table3x3_overlay_integratedAlpha_pTminCompare.png");

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
              const int nRows = 3;

              const int n10 = hPt10->GetXaxis()->GetNbins();
              const int n5  = hPt5->GetXaxis()->GetNbins();
              const int n3  = hPt3->GetXaxis()->GetNbins();
              const int nPt = std::min(n10, std::min(n5, n3));

              const std::string outPerPtDir = JoinPath(outRecoDir, "perPtBin_pTminCompare");
              EnsureDir(outPerPtDir);

              TCanvas c(TString::Format("c_tbl_pTminCompare_%s", rKey.c_str()).Data(),
                          TString::Format("c_tbl_pTminCompare_%s", rKey.c_str()).Data(),
                          1500, 1200);
              c.Divide(nCols, nRows, 0.001, 0.001);

              std::vector<TH1*> keep;
              keep.reserve(3 * 9);

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
              const int nPads   = std::min(9, nAvail);

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

                  const int b1 = (h->GetXaxis() ? h->GetXaxis()->GetFirst() : 1);
                  const int b2 = (h->GetXaxis() ? h->GetXaxis()->GetLast()  : h->GetNbinsX());

                  for (int ib = b1; ib <= b2; ++ib)
                  {
                    const double v = h->GetBinContent(ib) + h->GetBinError(ib);
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
                      if (ymaxRaw > 0.0) hSamRaw->SetMaximum(ymaxRaw * 1.10);
                      hSamRaw->SetMinimum(0.0);

                      const double yLineMax = (ymaxRaw > 0.0 ? ymaxRaw * 1.10 : 1.0);

                      TCanvas cRaw("c_SamVsJustin_raw","c_SamVsJustin_raw",900,700);
                      ApplyCanvasMargins1D(cRaw);

                        hSamRaw->Draw("E1");
                        hJustinRaw->Draw("E1 same");

                        // Cut floor lines (xJ,min bounds): abs and full
                        TLine* lAbs  = nullptr;
                        TLine* lFull = nullptr;
                        TLegend* legCuts = nullptr;

                        const double xAbs  = jetMinPtGeV / ptHi;
                        const double xFull = jetMinPtGeV / ptLo;

                        lAbs  = new TLine(xAbs,  0.0, xAbs,  yLineMax);
                        lFull = new TLine(xFull, 0.0, xFull, yLineMax);

                        lAbs->SetLineColor(kGreen+2);
                        lAbs->SetLineWidth(2);
                        lAbs->SetLineStyle(2);

                        lFull->SetLineColor(kOrange+7);
                        lFull->SetLineWidth(2);
                        lFull->SetLineStyle(2);

                        lAbs->Draw("same");
                        lFull->Draw("same");

                        // Title (centered)
                        TLatex tTitle;
                        tTitle.SetNDC(true);
                        tTitle.SetTextFont(42);
                        tTitle.SetTextAlign(22);
                        tTitle.SetTextSize(0.040);
                        tTitle.DrawLatex(0.50, 0.95, "Photon 10 + 20 GeV #gamma+Jet MC");

                      // Smaller legend, pushed to top-right (moved down to clear the centered title)
                      TLegend leg(0.52, 0.80, 0.86, 0.92);
                      leg.SetTextFont(42);
                      leg.SetTextSize(0.028);
                      leg.SetFillStyle(0);
                      leg.SetBorderSize(0);
                      leg.AddEntry(hSamRaw,    TString::Format("Sam's RECO (R = %.1f)", R).Data(), "ep");
                      leg.AddEntry(hJustinRaw, TString::Format("Justin's RECO w/ matched cuts (R = %.1f)", R).Data(), "ep");
                      leg.Draw();

                      // Cut-line legend under the main legend
                      legCuts = new TLegend(0.52, 0.70, 0.86, 0.80);
                      legCuts->SetTextFont(42);
                      legCuts->SetTextSize(0.028);
                      legCuts->SetFillStyle(0);
                      legCuts->SetBorderSize(0);
                      legCuts->AddEntry(lAbs,  "x_{J,min}^{abs} = p_{T}^{jet,min}/p_{T,max}^{#gamma}", "l");
                      legCuts->AddEntry(lFull, "x_{J,min}^{full} = p_{T}^{jet,min}/p_{T,min}^{#gamma}", "l");
                      legCuts->Draw();

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

                  double ymax = 0.0;
                  ymax = std::max(ymax, MaxBinContent(hSam));
                  ymax = std::max(ymax, MaxBinContent(hJustin));
                  if (ymax > 0.0) hSam->SetMaximum(ymax * 1.10);
                  hSam->SetMinimum(0.0);

                  const double yLineMax = (ymax > 0.0 ? ymax * 1.10 : 1.0);

                  TCanvas c("c_SamVsJustin","c_SamVsJustin",900,700);
                  ApplyCanvasMargins1D(c);

                    hSam->Draw("E1");
                  hJustin->Draw("E1 same");

                  // Cut floor lines (xJ,min bounds): abs and full
                  TLine* lAbs  = nullptr;
                  TLine* lFull = nullptr;
                  TLegend* legCuts = nullptr;

                  const double xAbs  = jetMinPtGeV / ptHi;
                  const double xFull = jetMinPtGeV / ptLo;

                  lAbs  = new TLine(xAbs,  0.0, xAbs,  yLineMax);
                  lFull = new TLine(xFull, 0.0, xFull, yLineMax);

                  lAbs->SetLineColor(kGreen+2);
                  lAbs->SetLineWidth(2);
                  lAbs->SetLineStyle(2);

                  lFull->SetLineColor(kOrange+7);
                  lFull->SetLineWidth(2);
                  lFull->SetLineStyle(2);

                  lAbs->Draw("same");
                  lFull->Draw("same");

                  // Title (centered)
                  TLatex tTitle;
                  tTitle.SetNDC(true);
                  tTitle.SetTextFont(42);
                  tTitle.SetTextAlign(22);
                  tTitle.SetTextSize(0.040);
                  tTitle.DrawLatex(0.50, 0.95, "Photon 10 + 20 GeV #gamma+Jet MC");

                  // Smaller legend, pushed to top-right (moved down to clear the centered title)
                  TLegend leg(0.52, 0.80, 0.86, 0.92);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.028);
                  leg.SetFillStyle(0);
                  leg.SetBorderSize(0);
                  leg.AddEntry(hSam,    TString::Format("Sam's RECO (R = %.1f)", R).Data(), "ep");
                  leg.AddEntry(hJustin, TString::Format("Justin's RECO w/ matched cuts (R = %.1f)", R).Data(), "ep");
                  leg.Draw();

                  // Cut-line legend under the main legend
                  legCuts = new TLegend(0.52, 0.70, 0.86, 0.80);
                  legCuts->SetTextFont(42);
                  legCuts->SetTextSize(0.028);
                  legCuts->SetFillStyle(0);
                  legCuts->SetBorderSize(0);
                  legCuts->AddEntry(lAbs,  "x_{J,min}^{abs} = p_{T}^{jet,min}/p_{T,max}^{#gamma}", "l");
                  legCuts->AddEntry(lFull, "x_{J,min}^{full} = p_{T}^{jet,min}/p_{T,min}^{#gamma}", "l");
                  legCuts->Draw();

                  // Info block: middle RHS
                  TLatex t;
                  t.SetNDC(true);
                  t.SetTextFont(42);
                  t.SetTextSize(0.032);
                  t.SetTextAlign(12);
                  t.DrawLatex(0.65, 0.62, TString::Format("p_{T}^{#gamma}: %.0f-%.0f GeV", ptLo, ptHi).Data());
                  t.DrawLatex(0.65, 0.57, TString::Format("Sam p_{T}^{#gamma} used: %s GeV", samPtLabel.c_str()).Data());
                  t.DrawLatex(0.65, 0.52, TString::Format("p_{T}^{jet,min} = %.0f GeV", jetMinPtGeV).Data());
                  t.DrawLatex(0.65, 0.47, TString::Format("Back-to-back: %s", bbLabel.c_str()).Data());

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

      #include "AnalyzeRecoilJets_RunJES3QA.cpp"


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

                    TAxis* ax = H.hRecoDphi->GetXaxis();
                    const int nPtUse = (ax ? ax->GetNbins() : 0);

                    for (int xbin = 1; xbin <= nPtUse; ++xbin)
                    {
                      const double axlo = ax->GetBinLowEdge(xbin);
                      const double axhi = ax->GetBinUpEdge(xbin);
                      const int lo = (int)std::lround(axlo);
                      const int hi = (int)std::lround(axhi);
                      const string folder = TString::Format("pT_%d_%d", lo, hi).Data();

                      TH1D* p = ProjY_AtXbin(H.hRecoDphi, xbin,
                        TString::Format("p_absDphi_reco_%s_%s_%s", ds.label.c_str(), rKey.c_str(), folder.c_str()).Data()
                      );
                      if (!p) continue;

                      const string perDir = JoinPath(dphiRecoDir, folder);
                      EnsureDir(perDir);

                      vector<string> lines = {
                        "RECO inclusive |#Delta#phi(#gamma,jet)|",
                        TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                        TString::Format("p_{T}^{#gamma}: %d-%d GeV", lo, hi).Data()
                      };

                      DrawAndSaveTH1_Common(ds, p,
                        JoinPath(perDir, "absDphi_counts.png"),
                        "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1"
                      );

                      TH1* ps = CloneTH1(p,
                        TString::Format("p_absDphi_reco_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), xbin).Data()
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

                    TAxis* ax = H.hTruthDphi->GetXaxis();
                    const int nPtUse = (ax ? ax->GetNbins() : 0);

                    for (int xbin = 1; xbin <= nPtUse; ++xbin)
                    {
                      const double axlo = ax->GetBinLowEdge(xbin);
                      const double axhi = ax->GetBinUpEdge(xbin);
                      const int lo = (int)std::lround(axlo);
                      const int hi = (int)std::lround(axhi);
                      const string folder = TString::Format("pT_%d_%d", lo, hi).Data();

                      TH1D* p = ProjY_AtXbin(H.hTruthDphi, xbin,
                        TString::Format("p_absDphi_truth_%s_%s_%s", ds.label.c_str(), rKey.c_str(), folder.c_str()).Data()
                      );
                      if (!p) continue;

                      const string perDir = JoinPath(dphiTruthDir, folder);
                      EnsureDir(perDir);

                      vector<string> lines = {
                        "TRUTH inclusive |#Delta#phi(#gamma,jet)|",
                        TString::Format("rKey=%s (R=%.1f)", rKey.c_str(), R).Data(),
                        TString::Format("p_{T}^{#gamma}: %d-%d GeV", lo, hi).Data()
                      };

                      DrawAndSaveTH1_Common(ds, p,
                        JoinPath(perDir, "absDphi_counts.png"),
                        "|#Delta#phi| [rad]", "Counts", lines, false, false, 0.0, "E1"
                      );

                      TH1* ps = CloneTH1(p,
                        TString::Format("p_absDphi_truth_shape_%s_%s_%d", ds.label.c_str(), rKey.c_str(), xbin).Data()
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
                const int nPtCanon = nPtReco;

                // Map canonical pT bin [lo,hi) to the unfolding axis bin by bin-center lookup.
                // NOTE: unfolding pT binning is allowed to be coarser (e.g. 10-15 as an "underflow-catch"
                // for canonical 13-15). We proceed as long as the canonical range is CONTAINED in the axis bin.
                auto MapCanonicalToAxisBin =
                  [&](TAxis* ax, const PtBin& b, const string& who)->int
                {
                  if (!ax) return -1;

                  const double mid  = 0.5 * (b.lo + b.hi);
                  const int    xbin = ax->FindBin(mid);

                  if (xbin < 1 || xbin > ax->GetNbins()) return -1;

                  const double axlo = ax->GetBinLowEdge(xbin);
                  const double axhi = ax->GetBinUpEdge(xbin);

                  const bool exact =
                    (std::fabs(axlo - (double)b.lo) <= 1e-6 && std::fabs(axhi - (double)b.hi) <= 1e-6);

                  const bool contained =
                    (axlo <= (double)b.lo + 1e-6 && axhi >= (double)b.hi - 1e-6);

                  if (!exact)
                  {
                    if (!contained)
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Unfolding pT mapping mismatch (" << who << "): requested "
                           << b.lo << "-" << b.hi
                           << " but axis bin " << xbin << " is "
                           << std::fixed << std::setprecision(0) << axlo << "-" << axhi
                           << " (NOT CONTAINED) → skipping this pT bin"
                           << ANSI_RESET << "\n";
                      return -1;
                    }

                    cout << ANSI_BOLD_YEL
                         << "[WARN] Unfolding pT mapping mismatch (" << who << "): requested "
                         << b.lo << "-" << b.hi
                         << " but axis bin " << xbin << " is "
                         << std::fixed << std::setprecision(0) << axlo << "-" << axhi
                         << " → proceeding using axis bin " << xbin
                         << ANSI_RESET << "\n";
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

                // Per unfolding pT bin (encoded on the unfolding axis): terminal line + text summary + overlay truth vs reco xJ shapes
                for (int xReco = 1; xReco <= nPtReco; ++xReco)
                {
                  const double axlo = H.hReco->GetXaxis()->GetBinLowEdge(xReco);
                  const double axhi = H.hReco->GetXaxis()->GetBinUpEdge(xReco);

                  PtBin b;
                  b.lo = (int)std::lround(axlo);
                  b.hi = (int)std::lround(axhi);

                  const int xTruth = MapCanonicalToAxisBin(H.hTruth->GetXaxis(), b, "truth");
                  if (xTruth < 1) continue;

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

                  // Truth-vs-reco xJ shape overlay per unfolding pT bin (aligned by edges)
                  {
                    TH1D* pxTruth = H.hTruth->ProjectionY(
                      TString::Format("pxTruth_%s_%d", rKey.c_str(), xReco).Data(), xTruth, xTruth
                    );
                    TH1D* pxReco  = H.hReco ->ProjectionY(
                      TString::Format("pxReco_%s_%d", rKey.c_str(), xReco).Data(),  xReco,  xReco
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
                        TString::Format("c_unf_ov_%s_%d", rKey.c_str(), xReco).Data(),
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

                      SaveCanvas(c, JoinPath(rOut, TString::Format("overlay_truth_vs_reco_xJ_shape_pTbin%d.png", xReco).Data()));

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
  
      #include "AnalyzeRecoilJets_RooUnfoldPipeline.cpp"
      #include "AnalyzeRecoilJets_ppAuAuComparisons.cpp"
      // =============================================================================
      // NEW: SIM+DATA PP — PPG12 SS template tables (Data vs Signal MC vs Background MC)
      //
      // Output:
      //   <kOutPPAuAuBase>/noIsoRequired/<pTbin>/table1x5_PP_SS_<tag>_DataSigBkg.png
      // where <tag> ∈ {pre, tight, nonTight}
      //
      // Notes:
      //   - Forces the SIM source for *_sig/*_bkg SS templates to the DEFAULT merged
      //     photonJet10+20 file keyed by kDefaultSimSampleKey (jetMinPt5_7piOver8).
      //   - Runs even when isPPdataAndAUAU == false (intended for isSimAndDataPP mode).
      // =============================================================================
      void RunPPG12SSTables_DataSigBkg_PP(const Dataset& dsPP, const Dataset& dsSIM)
      {
        cout << ANSI_BOLD_CYN
             << "\n[EXTRA] PPG12 SS template tables (SIM+DATA PP): Data vs Signal MC vs Background MC\n"
             << ANSI_RESET;

        if (!dsPP.topDir || !dsSIM.topDir)
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] PPG12 SS template tables skipped: PP or SIM dataset topDir is null."
               << ANSI_RESET << "\n";
          return;
        }

        const auto& ptBinsLocal = PtBins();
        if (ptBinsLocal.empty())
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] PPG12 SS template tables skipped: PtBins() is empty."
               << ANSI_RESET << "\n";
          return;
        }

        // Force the SIM input used for the overlays to be the DEFAULT merged photonJet10+20 file
        // keyed by kDefaultSimSampleKey (jetMinPt5_7piOver8).
        TFile* fSimSS = nullptr;
        TDirectory* simTopSS = dsSIM.topDir;

        if (bothPhoton10and20sim)
        {
            const string simMerged = MergedSimPath("photonJet10and20merged_SIM", "RecoilJets_photonjet10plus20_MERGED.root");

            cout << ANSI_DIM
                 << "  [SS templates] CfgTag()                      = " << CfgTag() << "\n"
                 << "  [SS templates] MergedSimPath(10+20)           = " << simMerged << "\n"
                 << ANSI_RESET;

          if (!simMerged.empty())
          {
            fSimSS = TFile::Open(simMerged.c_str(), "READ");
            if (fSimSS && !fSimSS->IsZombie())
            {
              TDirectory* d = fSimSS->GetDirectory(kDirSIM.c_str());
              if (d) simTopSS = d;
              else
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] SS templates: missing topDir '" << kDirSIM
                     << "' in merged SIM file: " << simMerged
                     << ANSI_RESET << "\n";
              }
            }
            else
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] SS templates: cannot open merged SIM file: " << simMerged
                   << ANSI_RESET << "\n";
              if (fSimSS) { fSimSS->Close(); delete fSimSS; }
              fSimSS = nullptr;
            }
          }
        }

        const string outBase = JoinPath(OutputPPAuAu(), "noIsoRequired");
        EnsureDir(outBase);

        struct VarDef { std::string var; std::string label; };
        const std::vector<VarDef> vars =
        {
          {"weta",   "w_{#eta}"},
          {"wphi",   "w_{#phi}"},
          {"e11e33", "E_{11}/E_{33}"},
          {"et1",    "et1"},
          {"e32e35", "E_{32}/E_{35}"}
        };

        auto LabelForVar = [&](const std::string& v) -> std::string
        {
          for (const auto& vd : vars) { if (vd.var == v) return vd.label; }
          return v;
        };

        const std::vector<std::string> ppg12Tags = {"pre", "tight", "nonTight"};

        for (int ipt = 0; ipt < (int)ptBinsLocal.size(); ++ipt)
        {
          const PtBin& pb = ptBinsLocal[ipt];
          const string outPt = JoinPath(outBase, pb.folder);
          EnsureDir(outPt);

          for (const auto& tag : ppg12Tags)
          {
            TCanvas cPP(
              TString::Format("c_pp_ss_%s_%s", tag.c_str(), pb.folder.c_str()).Data(),
              "c_pp_ss", 2600, 750
            );
            cPP.Divide(5, 1, 0.001, 0.001);

            std::vector<TH1*> keepAlive;
            keepAlive.reserve(vars.size() * 3);

            std::vector<TLegend*> keepLeg;
            keepLeg.reserve(vars.size());

            bool anyPad = false;

            for (int iv = 0; iv < (int)vars.size(); ++iv)
            {
              const std::string& var = vars[iv].var;
              const std::string  vlabel = LabelForVar(var);

              cPP.cd(iv + 1);
              gPad->SetLeftMargin(0.14);
              gPad->SetRightMargin(0.05);
              gPad->SetBottomMargin(0.14);
              gPad->SetTopMargin(0.18);
              gPad->SetLogy(false);

              const string hDataName = string("h_ss_") + var + string("_") + tag + pb.suffix;
              const string hSigName  = string("h_ss_") + var + string("_") + tag + string("_sig") + pb.suffix;
              const string hBkgName  = string("h_ss_") + var + string("_") + tag + string("_bkg") + pb.suffix;

              TH1* rawData = GetTH1FromTopDir(dsPP.topDir, hDataName);
              TH1* rawSig  = GetTH1FromTopDir(simTopSS,  hSigName);
              TH1* rawBkg  = GetTH1FromTopDir(simTopSS,  hBkgName);

              if (!rawData && !rawSig && !rawBkg)
              {
                DrawMissingPad(TString::Format("%s, %s, %s", var.c_str(), tag.c_str(), pb.folder.c_str()).Data());
                continue;
              }

              anyPad = true;

                TH1* hData = nullptr;
                TH1* hSig  = nullptr;
                TH1* hBkg  = nullptr;

                if (rawData)
                {
                  hData = CloneNormalizeStyle(rawData,
                    TString::Format("ss_%s_%s_%s_data", tag.c_str(), var.c_str(), pb.folder.c_str()).Data(),
                    kBlack, 20);

                  if (hData)
                  {
                    hData->SetLineWidth(2);
                    hData->SetLineColor(kBlack);
                    hData->SetMarkerColor(kBlack);
                    hData->SetMarkerStyle(20);
                    hData->SetMarkerSize(1.00);
                    hData->SetFillStyle(0);
                  }
                }

                if (rawSig)
                {
                  hSig = CloneNormalizeStyle(rawSig,
                    TString::Format("ss_%s_%s_%s_sig", tag.c_str(), var.c_str(), pb.folder.c_str()).Data(),
                    kRed + 1, 24);

                  if (hSig)
                  {
                    hSig->SetLineWidth(2);
                    hSig->SetLineColor(kRed + 1);
                    hSig->SetMarkerColor(kRed + 1);
                    hSig->SetMarkerStyle(1);
                    hSig->SetMarkerSize(0.0);
                    hSig->SetFillStyle(0);
                  }
                }

                if (rawBkg)
                {
                  hBkg = CloneNormalizeStyle(rawBkg,
                    TString::Format("ss_%s_%s_%s_bkg", tag.c_str(), var.c_str(), pb.folder.c_str()).Data(),
                    kBlue + 1, 25);

                  if (hBkg)
                  {
                    hBkg->SetLineWidth(2);
                    hBkg->SetLineColor(kBlue + 1);
                    hBkg->SetMarkerColor(kBlue + 1);
                    hBkg->SetMarkerStyle(1);
                    hBkg->SetMarkerSize(0.0);
                    hBkg->SetFillStyle(0);
                  }
                }

                TH1* hFirst = (hData ? hData : (hSig ? hSig : hBkg));
                if (!hFirst)
                {
                  DrawMissingPad(TString::Format("%s, %s, %s", var.c_str(), tag.c_str(), pb.folder.c_str()).Data());
                  continue;
                }

                TH1* hFrame = (hSig ? hSig : (hBkg ? hBkg : hData));

                hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                hFrame->GetYaxis()->SetTitle("Unit Normalized");

                double yMax = 0.0;
                if (hData)
                {
                  for (int ib = 1; ib <= hData->GetNbinsX(); ++ib)
                  {
                    yMax = std::max(yMax, (double)(hData->GetBinContent(ib) + hData->GetBinError(ib)));
                  }
                }
                if (hSig)  yMax = std::max(yMax, (double)hSig->GetMaximum());
                if (hBkg)  yMax = std::max(yMax, (double)hBkg->GetMaximum());

                const double yScale = ((var == "weta" || var == "wphi") ? 1.35 : 1.35);

                hFrame->SetMinimum(0.0);
                hFrame->SetMaximum((yMax > 0.0) ? (yMax * yScale) : 1.0);

                if (hSig)
                {
                  hSig->Draw("HIST");
                }
                else if (hBkg)
                {
                  hBkg->Draw("HIST");
                }
                else
                {
                  hData->Draw("E1");
                }

                if (hBkg && hBkg != hSig)  hBkg->Draw("HIST same");
                if (hData)                 hData->Draw("E1 same");

                // Legend
                {
                  const bool isW = (var == "weta" || var == "wphi");
                  TLegend* leg = (isW ? new TLegend(0.55, 0.61, 0.93, 0.8) : new TLegend(0.16, 0.61, 0.54, 0.8));
                  leg->SetBorderSize(0);
                  leg->SetFillStyle(0);
                  leg->SetTextFont(42);
                  leg->SetTextSize(0.038);

                  if (hData) leg->AddEntry(hData, "Data", "ep");
                  if (hSig)  leg->AddEntry(hSig,  "Signal MC", "l");
                  if (hBkg)  leg->AddEntry(hBkg,  "Background MC", "l");

                  leg->Draw();
                  keepLeg.push_back(leg);
                }

                // Pad header
                {
                  std::string tagLabel = "Preselection";
                  if (tag == "tight") tagLabel = "Tight";
                  else if (tag == "nonTight") tagLabel = "Non-tight";

                  TLatex th;
                  th.SetNDC(true);
                  th.SetTextFont(42);
                  th.SetTextAlign(22);
                  th.SetTextSize(0.050);
                  th.DrawLatex(0.50, 0.91,
                    TString::Format("%s, %s, p_{T}^{#gamma}: %d-%d GeV",
                      vlabel.c_str(), tagLabel.c_str(), pb.lo, pb.hi).Data());
                }

                // SS cut label + cut lines (fixed-cut vars)
                {
                  TLatex tcut;
                  tcut.SetNDC(true);
                  tcut.SetTextFont(42);
                  tcut.SetTextAlign(13);
                  tcut.SetTextSize(0.040);

                  bool drawCuts = false;
                  bool drawSingleCut = false;
                  double cutLo = 0.0;
                  double cutHi = 0.0;
                  std::string cutText;

                  if (var == "e11e33")
                  {
                    cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                    drawCuts = true;
                    cutLo = 0.4;
                    cutHi = 0.98;
                  }
                  else if (var == "e32e35")
                  {
                    cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                    drawCuts = true;
                    cutLo = 0.92;
                    cutHi = 1.0;
                  }
                  else if (var == "et1")
                  {
                    cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                    drawCuts = true;
                    cutLo = 0.9;
                    cutHi = 1.0;
                  }
                  else if (var == "weta")
                  {
                    cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    drawSingleCut = true;
                    const double ptCenter = 0.5 * (pb.lo + pb.hi);
                    cutHi = 0.15 + 0.006 * ptCenter;
                  }
                  else if (var == "wphi")
                  {
                    cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                    drawSingleCut = true;
                    const double ptCenter = 0.5 * (pb.lo + pb.hi);
                    cutHi = 0.15 + 0.006 * ptCenter;
                  }

                  if (!cutText.empty())
                  {
                    tcut.DrawLatex(0.16, 0.86, cutText.c_str());
                  }

                  if (drawCuts || drawSingleCut)
                  {
                    gPad->Update();
                    const double yMin = gPad->GetUymin();
                    const double yMaxPad = gPad->GetUymax();

                    if (drawCuts)
                    {
                      TLine* l1 = new TLine(cutLo, yMin, cutLo, yMaxPad);
                      l1->SetLineColor(kBlack);
                      l1->SetLineWidth(2);
                      l1->SetLineStyle(2);
                      l1->Draw("same");

                      TLine* l2 = new TLine(cutHi, yMin, cutHi, yMaxPad);
                      l2->SetLineColor(kBlack);
                      l2->SetLineWidth(2);
                      l2->SetLineStyle(2);
                      l2->Draw("same");
                    }

                    if (drawSingleCut)
                    {
                      TLine* l1 = new TLine(cutHi, yMin, cutHi, yMaxPad);
                      l1->SetLineColor(kBlack);
                      l1->SetLineWidth(2);
                      l1->SetLineStyle(2);
                      l1->Draw("same");
                    }
                  }

                  gPad->RedrawAxis();
                }

              if (hData) keepAlive.push_back(hData);
              if (hSig)  keepAlive.push_back(hSig);
              if (hBkg)  keepAlive.push_back(hBkg);
            }

            if (anyPad)
            {
              SaveCanvas(cPP, JoinPath(outPt,
                TString::Format("table1x5_PP_SS_%s_DataSigBkg.png", tag.c_str()).Data()));
            }

            for (TLegend* l : keepLeg) delete l;
            keepLeg.clear();

            for (TH1* h : keepAlive) delete h;
            keepAlive.clear();
          }
        }

        if (fSimSS)
        {
          fSimSS->Close();
          delete fSimSS;
        }
      }

      // =============================================================================
      // High-level runner
      // =============================================================================
      // NOTE: Run-mode validation now lives in AnalyzeRecoilJets.h (ValidateRunConfig()).
      // ExactlyOneModeSet() is provided by the header (ARJ scope) and no longer
      // needs a local wrapper here.

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
        if (gPPOutBaseSubdir.empty()) return OutputPP();
        return JoinPath(OutputPP(), gPPOutBaseSubdir);
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

        const bool doAuAu = (mode == RunMode::kAuAuOnly);

        const SimSample ss = CurrentSimSample();

        if (doSim)
        {
            if (ss == SimSample::kSimEmbedded)
            {
              // Embedded SIM: per-centrality datasets (same file, filtered by centSuffix)
              // For embedded SIM, derive the cfg tag from the actual ROOT file name so the
              // filename itself is the source of truth for vz / cone / iso / UE variant.
              //
              // IMPORTANT:
              // The current embedded photon20 files were only analyzed for:
              //   0-10, 10-20, 20-40, 40-60, 60-80
              // and do NOT contain 80-100 histograms. So suppress 80-100 here to avoid
              // creating empty output folders like .../80_100.
              vector<CentBin> centBins;
              for (const auto& cb : CentBins())
              {
                if (cb.lo == 80 && cb.hi == 100) continue;
                centBins.push_back(cb);
              }

                const string simInPath  = SimInputPathForSample(ss);
                const string simCfgTag  = CfgTagWithUE();
                const string simOutBase = OutputSimEmbedded();

                if (centBins.empty())
                {
                  Dataset ds;
                  ds.label      = "SIM_EMBEDDED_" + simCfgTag;
                  ds.isSim      = true;
                  ds.trigger    = "";
                  ds.topDirName = kDirSIM;
                  ds.inFilePath = simInPath;
                  ds.outBase    = simOutBase;
                  datasets.push_back(std::move(ds));
              }
              else
              {
                for (const auto& cb : centBins)
                {
                  Dataset ds;
                  ds.label      = TString::Format("SIM_EMBEDDED_%s_%d_%d", simCfgTag.c_str(), cb.lo, cb.hi).Data();
                  ds.isSim      = true;
                  ds.trigger    = "";
                  ds.topDirName = kDirSIM;
                  ds.inFilePath = simInPath;
                  ds.centFolder = cb.folder;
                  ds.centSuffix = cb.suffix;
                  ds.centLabel  = TString::Format("Centrality: %d-%d%%", cb.lo, cb.hi).Data();
                  ds.outBase    = JoinPath(simOutBase, cb.folder);
                  datasets.push_back(std::move(ds));
                }
              }
            }
          else
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
        }

        if (doPP)
        {
          Dataset ds;
          ds.label      = "DATA_PP";
          ds.isSim      = false;
          ds.trigger    = kTriggerPP;
          ds.topDirName = kTriggerPP;
          ds.inFilePath = InputPP(isRun25pp);

            // NOTE:
            // - Legacy single-run behavior: PP output goes to OutputPP()
            // - Multi-run SIM+DATA behavior: PP output goes to OutputPP()/with_<simSampleLabel>
            //   (so you can loop over multiple SIM selections without clobbering PP outputs)
            //
            // DATA outputs are organized per-trigger:
            //   <PP base>/<trigger>/{baselineData,insituCalib,unfolding,...}
            if (mode == RunMode::kPPDataOnly)
            {
              ds.outBase = JoinPath(OutputPP(), ds.trigger);
            }
            else
            {
              ds.outBase = JoinPath(PPOutBaseForThisRun(), ds.trigger);
            }

          datasets.push_back(std::move(ds));
        }

        if (doAuAu)
        {
          const auto& centBins = CentBins();

          if (centBins.empty())
          {
            Dataset ds;
            ds.label      = "DATA_AUAU";
            ds.isSim      = false;
            ds.trigger    = kTriggerAuAuGold;
            ds.topDirName = kTriggerAuAuGold;
            ds.inFilePath = InputAuAu();

            // Fallback: no YAML centrality bins available, keep legacy AuAu output path.
            ds.outBase = JoinPath(OutputAuAu(), ds.trigger);

            datasets.push_back(std::move(ds));
          }
          else
          {
            for (const auto& cb : centBins)
            {
              Dataset ds;
              ds.label      = TString::Format("DATA_AUAU_%d_%d", cb.lo, cb.hi).Data();
              ds.isSim      = false;
              ds.trigger    = kTriggerAuAuGold;
              ds.topDirName = kTriggerAuAuGold;
              ds.inFilePath = InputAuAu();

              ds.centFolder = cb.folder;
              ds.centSuffix = cb.suffix;
              ds.centLabel  = TString::Format("Centrality: %d-%d%%", cb.lo, cb.hi).Data();

              // AuAu-only outputs are organized per-trigger, per-centrality:
              //   <AuAu base>/<trigger>/<centFolder>/{baselineData,insituCalib,unfolding,...}
              ds.outBase = JoinPath(JoinPath(OutputAuAu(), ds.trigger), ds.centFolder);

              datasets.push_back(std::move(ds));
            }
          }
        }

        return datasets;
    }

    inline bool MaybeBuildMergedSIM(RunMode mode)
    {
        // Only relevant when a SIM-including mode is running AND a merged SIM sample was selected.
        if (mode == RunMode::kPPDataOnly || mode == RunMode::kAuAuOnly) return true;

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
              const auto& cfg = DefaultSim10and20Config();

              const string outMerged =
                  MergedSIMOut_10and20_ForKey(cfg.key);

              cout << ANSI_BOLD_CYN
                   << "\n[MERGE SIM] Rebuilding ONLY the default SIM10+20 merged file.\n"
                   << "            cfg.key = " << cfg.key << "\n"
                   << "            in10    = " << cfg.photon10 << "\n"
                   << "            in20    = " << cfg.photon20 << "\n"
                   << "            out     = " << outMerged << "\n"
                   << ANSI_RESET;

              ok = BuildMergedSIMFile_PhotonSlices(
                {cfg.photon10, cfg.photon20},
                {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                outMerged,
                kDirSIM,
                {"photonJet10", "photonJet20"}
              );
          }
      }
      else if (ss == SimSample::kPhotonJet5And10And20Merged)
      {
          if (!doRemergePhoton5and10and20sim)
          {
              const string outMerged = kMergedSIMOut_5and10and20;

              cout << ANSI_BOLD_CYN
                   << "\n[MERGE SIM] doRemergePhoton5and10and20sim=false -> skipping SIM5+10+20 rebuild step.\n"
                   << "            Using existing merged output: " << outMerged << "\n"
                   << ANSI_RESET;

              if (gSystem->AccessPathName(outMerged.c_str()))
              {
                cout << ANSI_BOLD_RED
                     << "[MERGE SIM][FATAL] Merged SIM5+10+20 file not found, but doRemergePhoton5and10and20sim=false:\n"
                     << "  " << outMerged << "\n"
                     << ANSI_RESET;
                ok = false;
              }
          }
          else
          {
              const auto& cfg = DefaultSim10and20Config();

              const string outMerged =
                  kMergedSIMOut_5and10and20;

              cout << ANSI_BOLD_CYN
                   << "\n[MERGE SIM] Rebuilding ONLY the default SIM5+10+20 merged file.\n"
                   << "            cfg.key = " << cfg.key << "\n"
                   << "            in5     = " << kInSIM5 << "\n"
                   << "            in10    = " << cfg.photon10 << "\n"
                   << "            in20    = " << cfg.photon20 << "\n"
                   << "            out     = " << outMerged << "\n"
                   << ANSI_RESET;

              ok = BuildMergedSIMFile_PhotonSlices(
                {kInSIM5, cfg.photon10, cfg.photon20},
                {kSigmaPhoton5_pb, kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                outMerged,
                kDirSIM,
                {"photonJet5", "photonJet10", "photonJet20"}
              );
          }
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
      static bool sDidLoadSPhenixStyle = false;
      if (!sDidLoadSPhenixStyle)
      {
          sDidLoadSPhenixStyle = true;
          if (!gROOT->GetStyle("sPHENIX"))
          {
            gROOT->LoadMacro("sPhenixStyle.C");
            gROOT->ProcessLine("SetsPhenixStyle();");
          }
          else
          {
            gROOT->SetStyle("sPHENIX");
            gROOT->ForceStyle();
          }
      }

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
           << "    ARJ_HAVE_ROOUNFOLD       = " << ARJ_HAVE_ROOUNFOLD << "\n"
           << "    isPPdataOnly            = " << (isPPdataOnly ? "true" : "false") << "\n"
           << "    isSimAndDataPP          = " << (isSimAndDataPP ? "true" : "false") << "\n"
           << "    isAuAuOnly              = " << (isAuAuOnly ? "true" : "false") << "\n"
           << "    isPhotonJet5            = " << (isPhotonJet5 ? "true" : "false") << "\n"
           << "    isPhotonJet10           = " << (isPhotonJet10 ? "true" : "false") << "\n"
           << "    isPhotonJet20           = " << (isPhotonJet20 ? "true" : "false") << "\n"
           << "    bothPhoton5and10sim     = " << (bothPhoton5and10sim ? "true" : "false") << "\n"
           << "    bothPhoton5and20sim     = " << (bothPhoton5and20sim ? "true" : "false") << "\n"
           << "    bothPhoton10and20sim    = " << (bothPhoton10and20sim ? "true" : "false") << "\n"
           << "    doRemergePhoton10and20sim = " << (doRemergePhoton10and20sim ? "true" : "false") << "\n"
           << "    allPhoton5and10and20sim = " << (allPhoton5and10and20sim ? "true" : "false") << "\n"
           << "    isSimEmbedded           = " << (isSimEmbedded ? "true" : "false") << "\n"
           << "      kEmbeddedCutTag       = " << kEmbeddedCutTag << "\n"
           << "      kEmbeddedUEVariant    = " << kEmbeddedUEVariant << "\n"
           << "      EmbeddedSimCfgTag()   = " << EmbeddedSimCfgTag() << "\n";

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
      if (mode == RunMode::kSimOnly || mode == RunMode::kSimAndDataPP)
      {
        cout << "  |  SIM sample: " << SimSampleLabel(ss)
             << "  |  SIM outBase: " << SimOutBaseForSample(ss);
      }
      cout << ANSI_RESET << "\n";

      // ---------------------------------------------------------------------------
      // Explicit default SIM config printout (so it's obvious what drives the merge + analysis)
      // ---------------------------------------------------------------------------
      if (mode == RunMode::kSimOnly || mode == RunMode::kSimAndDataPP)
      {
        const Sim10and20Config& cfgDef = DefaultSim10and20Config();

        cout << ANSI_DIM
             << "\n  [DEFAULT SIM CONFIG]\n"
             << "    DefaultSimSampleKey()           = " << DefaultSimSampleKey() << "\n"
             << "    cfg.key                         = " << cfgDef.key << "\n";

        if (ss == SimSample::kPhotonJet5And10And20Merged)
        {
          cout << "    cfg.photon5                     = " << cfgDef.photon5 << "\n";
        }

        cout << "    cfg.photon10                    = " << cfgDef.photon10 << "\n"
             << "    cfg.photon20                    = " << cfgDef.photon20 << "\n"
             << "    cfg.jetMinPt                    = " << cfgDef.jetMinPt << " GeV\n"
             << "    cfg.bbLabel                     = " << cfgDef.bbLabel << "\n";

        if (ss == SimSample::kPhotonJet5And10And20Merged)
        {
          cout << "    kInSIM5 (alias)                 = " << kInSIM5 << "\n";
        }

        cout << "    kInSIM10 (alias)                = " << kInSIM10 << "\n"
             << "    kInSIM20 (alias)                = " << kInSIM20 << "\n";

        if (ss == SimSample::kPhotonJet5And10And20Merged)
        {
          cout << "    kMergedSIMOut_5and10and20       = " << kMergedSIMOut_5and10and20 << "\n";
        }
        else
        {
          cout << "    MergedSIMOut_10and20_Default()  = " << MergedSIMOut_10and20_Default() << "\n";
        }

        cout << ANSI_RESET;

        if (ss == SimSample::kPhotonJet10And20Merged ||
            ss == SimSample::kPhotonJet5And10And20Merged)
        {
          cout << ANSI_DIM
               << "  [ACTIVE MERGED DATASET PATH]\n"
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
      // AuAu-only: Tabulate accepted events per centrality (and total) BEFORE any QA
      // ---------------------------------------------------------------------------
      if (isAuAuOnly)
      {
        cout << ANSI_BOLD_CYN << "\n[AuAuOnly] Accepted events summary (by centrality + total)\n" << ANSI_RESET;

        std::map<std::string, std::vector<std::pair<std::pair<int,int>, double>>> accByTrig;

        auto ParseCentLoHiFromSuffix = [&](const std::string& s, int& lo, int& hi)->bool
        {
          lo = -1; hi = -1;
          if (s.empty()) return false;

          int a = -1, b = -1;
          if (std::sscanf(s.c_str(), "_cent_%d_%d", &a, &b) == 2)
          {
            lo = a; hi = b;
            return true;
          }
          return false;
        };

        for (auto& ds : datasets)
        {
          if (ds.isSim) continue;
          if (ds.centSuffix.empty()) continue;

          int clo = -1, chi = -1;
          if (!ParseCentLoHiFromSuffix(ds.centSuffix, clo, chi)) continue;

          const double nEvt = ReadEventCount(ds);
          accByTrig[ds.trigger].push_back({{clo, chi}, nEvt});
        }

          for (auto& kv : accByTrig)
          {
            const std::string& trig = kv.first;
            auto& v = kv.second;

            if (v.empty()) continue;

            std::sort(v.begin(), v.end(),
                      [](const auto& a, const auto& b){ return a.first.first < b.first.first; });

            cout << ANSI_BOLD_YEL << "\n  Trigger: " << trig << ANSI_RESET << "\n";

            double total = 0.0;
            for (const auto& it : v)
            {
              const int clo = it.first.first;
              const int chi = it.first.second;
              const double nEvt = it.second;

              total += nEvt;

              cout << "    [DATA_AUAU_" << clo << "_" << chi << "] accepted events = "
                   << std::fixed << std::setprecision(0) << nEvt << "\n";
            }

            cout << "    [TOTAL] accepted events (sum over centrality) = "
                 << std::fixed << std::setprecision(0) << total << "\n";

            // ---------------------------------------------------------------------
            // Also save the centrality distribution histogram to the SAME trigger dir
            //   <kOutAuAuBase>/<trigger>/centrality_distribution.png
            // ---------------------------------------------------------------------
            Dataset* dsRef = nullptr;
            for (auto& ds : datasets)
            {
              if (ds.isSim) continue;
              if (ds.trigger != trig) continue;
              dsRef = &ds;
              break;
            }

            if (dsRef && dsRef->topDir)
            {
              TH1* hCent = dynamic_cast<TH1*>(dsRef->topDir->Get("h_centrality"));

              if (hCent)
              {
                const std::string outDir  = JoinPath(OutputAuAu(), trig);
                const std::string outPath = JoinPath(outDir, "centrality_distribution.png");
                EnsureDir(outDir);

                TCanvas c("cCentralityDist", "", 900, 700);
                c.SetTopMargin(0.12);
                c.SetBottomMargin(0.14);
                c.SetLeftMargin(0.13);
                c.SetRightMargin(0.05);

                hCent->SetTitle("");
                hCent->SetStats(0);
                hCent->GetXaxis()->SetTitle("Centrality [%]");
                hCent->GetYaxis()->SetTitle("Counts");
                hCent->GetYaxis()->SetTitleOffset(1.20);

                hCent->Draw("hist");

                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(23); // center, top
                t.SetTextSize(0.045);
                t.DrawLatex(0.50, 0.98, "Run25auau Centrality Distribution, MBD NS #geq 2 vtx < 150 cm");

                c.SaveAs(outPath.c_str());
              }
              else
              {
                cout << "    [WARN] h_centrality not found for trigger '" << trig << "'\n";
              }
            }
            else
            {
              cout << "    [WARN] Cannot save centrality plot (dataset not found/opened for trigger '" << trig << "')\n";
            }
          }
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

        cout << "  -> [pi0 QA] corrected vs no-asinh-correction overlays (DATA only)...\n";
        analysis::RunPi0QA(ds);
        cout << "     [OK] pi0 QA complete.\n";

        cout << "  -> [Section 2] Preselection failure table...\n";
        analysis::RunPreselectionFailureTable(ds);
        cout << "     [OK] Preselection failure table complete.\n";

        cout << "  -> [Section 3] General isolation QA...\n";
        analysis::RunIsolationQA(ds);
        cout << "     [OK] Isolation QA complete.\n";
      }

      // ---------------------------------------------------------------------------
      // AuAu-only: Accepted events vs centrality (one plot per trigger)
      //   Outputs to:  <kOutAuAuBase>/<trigger>/acceptedEvents_vs_centrality.png
      // ---------------------------------------------------------------------------
      if (isAuAuOnly)
      {
        std::map<std::string, std::vector<std::pair<std::pair<int,int>, double>>> accByTrig;

        auto ParseCentLoHiFromSuffix = [&](const std::string& s, int& lo, int& hi)->bool
        {
          lo = -1; hi = -1;
          if (s.empty()) return false;

          // Expect: "_cent_<lo>_<hi>"
          int a = -1, b = -1;
          if (std::sscanf(s.c_str(), "_cent_%d_%d", &a, &b) == 2)
          {
            lo = a; hi = b;
            return true;
          }
          return false;
        };

        for (auto& ds : datasets)
        {
          if (ds.isSim) continue;
          if (ds.centSuffix.empty()) continue;

          int clo = -1, chi = -1;
          if (!ParseCentLoHiFromSuffix(ds.centSuffix, clo, chi)) continue;

          const double nEvt = ReadEventCount(ds);
          accByTrig[ds.trigger].push_back({{clo, chi}, nEvt});
        }

        for (auto& kv : accByTrig)
        {
          const std::string& trig = kv.first;
          auto& v = kv.second;

          if (v.empty()) continue;

          std::sort(v.begin(), v.end(),
                    [](const auto& a, const auto& b){ return a.first.first < b.first.first; });

          const int n = (int)v.size();

          const std::string outDir  = JoinPath(OutputAuAu(), trig);
          const std::string outPath = JoinPath(outDir, "acceptedEvents_vs_centrality.png");
          EnsureDir(outDir);

          TH1D* h = new TH1D("hAcceptedEventsVsCent", "", n, 0.0, (double)n);
          h->SetStats(0);
          h->GetYaxis()->SetTitle("Accepted Events");
          h->GetXaxis()->SetTitle("Centrality [%]");
          h->GetXaxis()->SetLabelSize(0.045);
          h->GetYaxis()->SetTitleOffset(1.20);

          for (int i = 0; i < n; ++i)
          {
            const int clo = v[i].first.first;
            const int chi = v[i].first.second;
            const double nEvt = v[i].second;

            h->SetBinContent(i + 1, nEvt);
            h->GetXaxis()->SetBinLabel(i + 1, TString::Format("%d-%d", clo, chi).Data());
          }

          TCanvas* c = new TCanvas("cAcceptedEventsVsCent", "", 900, 650);
          c->SetTopMargin(0.12);
          c->SetBottomMargin(0.15);
          c->SetLeftMargin(0.14);
          c->SetRightMargin(0.05);

          h->Draw("hist");

          TLatex t;
          t.SetNDC(true);
          t.SetTextFont(42);
          t.SetTextAlign(23);   // center, top
          t.SetTextSize(0.045);
          t.DrawLatex(0.50, 0.98, "Run25auau Accepted Events, MBD NS #geq 2 vtx < 150 cm");

          c->SaveAs(outPath.c_str());

          delete c;
          delete h;
        }
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
      // [5I] RooUnfold pipeline (SIM+DATA PP only): unfold photons + (pTgamma,xJ) and produce per-photon xJ tables
      // ---------------------------------------------------------------------------
      {
        const bool simAndDataPPMergedForUnfold =
          (bothPhoton10and20sim || allPhoton5and10and20sim);

        cout << ANSI_BOLD_CYN
             << "\n[5I] RooUnfold gate check\n"
             << "  ARJ_HAVE_ROOUNFOLD        = " << ARJ_HAVE_ROOUNFOLD << "\n"
             << "  isSimAndDataPP            = " << (isSimAndDataPP ? "true" : "false") << "\n"
             << "  bothPhoton10and20sim      = " << (bothPhoton10and20sim ? "true" : "false") << "\n"
             << "  allPhoton5and10and20sim   = " << (allPhoton5and10and20sim ? "true" : "false") << "\n"
             << "  datasets.size()           = " << datasets.size() << "\n"
             << ANSI_RESET;

        for (auto& ds : datasets)
        {
          cout << "    - " << ds.label
               << "  isSim=" << (ds.isSim ? "true" : "false")
               << "  topDirName=" << ds.topDirName
               << "  inFilePath=" << ds.inFilePath
               << "  outBase=" << ds.outBase
               << "\n";
        }

          if (isSimAndDataPP && simAndDataPPMergedForUnfold &&
              (do_xJ_PPunfold || (!do_xJ_PPunfold && gApplyPurityCorrectionForUnfolding)))
          {
              Dataset* dsSIM = nullptr;
              Dataset* dsPP  = nullptr;

              for (auto& ds : datasets)
              {
                if (ds.isSim) dsSIM = &ds;
                else         dsPP  = &ds;
              }

              if (!dsSIM || !dsPP)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] RooUnfold pipeline requested, but SIM or DATA dataset is missing. Skipping."
                     << ANSI_RESET << "\n";
              }
              else
              {
                const bool runOnlyPurityCorrected = (!do_xJ_PPunfold && gApplyPurityCorrectionForUnfolding);

                if (!runOnlyPurityCorrected)
                {
                  cout << "  -> [5I] RooUnfold pipeline (SIM+DATA PP): non-purity-corrected unfold to particle level + per-photon x_{J} tables...\n";
                  gApplyPurityCorrectionForUnfolding = false;
                  analysis::RunRooUnfoldPipeline_SimAndDataPP(*dsPP, *dsSIM);
                  cout << "     [OK] nonPurityCorrected\n";
                }

                cout << "  -> [5I] RooUnfold pipeline (SIM+DATA PP): purity-corrected unfold to particle level + per-photon x_{J} tables...\n";
                gApplyPurityCorrectionForUnfolding = true;
                analysis::RunRooUnfoldPipeline_SimAndDataPP(*dsPP, *dsSIM);
                cout << "     [OK] purityCorrected\n";

                if (!runOnlyPurityCorrected)
                {
                  gApplyPurityCorrectionForUnfolding = false;
                  cout << "  -> [5I] purity-corrected vs non-purity-corrected per-photon x_{J} overlays...\n";
                  analysis::RunPurityCorrectedUncorrectedOverlayPP(*dsPP);
                  cout << "     [OK] purityCorrectedUncorrectedOverly\n";
                }
              }
            }
            else
            {
              cout << ANSI_BOLD_YEL
                   << "[5I] Skipping RooUnfold pipeline: requires isSimAndDataPP && (bothPhoton10and20sim || allPhoton5and10and20sim) and either do_xJ_PPunfold=true or (do_xJ_PPunfold=false with gApplyPurityCorrectionForUnfolding=true).\n"
                   << "     Current: isSimAndDataPP=" << (isSimAndDataPP ? "true" : "false")
                   << " bothPhoton10and20sim=" << (bothPhoton10and20sim ? "true" : "false")
                   << " allPhoton5and10and20sim=" << (allPhoton5and10and20sim ? "true" : "false")
                   << " do_xJ_PPunfold=" << (do_xJ_PPunfold ? "true" : "false")
                   << " gApplyPurityCorrectionForUnfolding=" << (gApplyPurityCorrectionForUnfolding ? "true" : "false")
                   << ANSI_RESET << "\n";
            }
      }

      // ---------------------------------------------------------------------------
      // [5J] PPG12 SS template tables (SIM+DATA PP): Data vs Signal MC vs Background MC
      // Output:
      //   <kOutPPAuAuBase>/noIsoRequired/<pTbin>/table1x5_PP_SS_<tag>_DataSigBkg.png
      // ---------------------------------------------------------------------------
      if (isSimAndDataPP)
      {
        Dataset* dsSIM = nullptr;
        Dataset* dsPP  = nullptr;

        for (auto& ds : datasets)
        {
          if (ds.isSim) dsSIM = &ds;
          else         dsPP  = &ds;
        }

        if (!dsSIM || !dsPP)
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] PPG12 SS template tables requested (isSimAndDataPP), but SIM or DATA dataset is missing. Skipping."
               << ANSI_RESET << "\n";
        }
        else
        {
          analysis::RunPPG12SSTables_DataSigBkg_PP(*dsPP, *dsSIM);
        }
      }

      // ---------------------------------------------------------------------------
      // OPTIONAL: PP vs Au+Au (gold-gold) photon-ID deliverables (requires both PP and AuAu)
      // ---------------------------------------------------------------------------
      if (isPPdataAndAUAU)
      {
          // Find the PP dataset: must match kTriggerPP (not just !isSim, which also matches AuAu)
          Dataset* dsPP = nullptr;
          for (auto& ds : datasets)
              {
                if (!ds.isSim && ds.trigger == kTriggerPP)
                {
                  dsPP = &ds;
                  break;
                }
          }

          // In kAuAuOnly mode the PP dataset is not in the datasets vector.
          // Open it on-the-fly so the PP vs AuAu overlays use the real PP data.
          Dataset ppFallback;
          TFile* ppFallbackFile = nullptr;
          if (!dsPP)
          {
                ppFallbackFile = TFile::Open(InputPP(isRun25pp).c_str(), "READ");
                if (ppFallbackFile && !ppFallbackFile->IsZombie())
                {
                  ppFallback.label      = "DATA_PP";
                  ppFallback.isSim      = false;
                  ppFallback.trigger    = kTriggerPP;
                  ppFallback.topDirName = kTriggerPP;
                  ppFallback.inFilePath = InputPP(isRun25pp);
                  ppFallback.outBase    = JoinPath(OutputPP(), kTriggerPP);
                  ppFallback.file       = ppFallbackFile;
                  ppFallback.topDir     = ppFallbackFile->GetDirectory(kTriggerPP.c_str());

                  if (ppFallback.topDir)
                  {
                    dsPP = &ppFallback;
                    cout << ANSI_BOLD_CYN
                         << "[INFO] Opened PP file on-the-fly for PP vs AuAu overlays: " << InputPP(isRun25pp)
                         << ANSI_RESET << "\n";
                  }
                  else
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] PP file opened but missing trigger dir '" << kTriggerPP
                         << "' in: " << InputPP(isRun25pp) << ". Skipping PP vs AuAu overlays."
                         << ANSI_RESET << "\n";
                  }
                }
                else
                {
                  cout << ANSI_BOLD_YEL
                       << "[WARN] isPPdataAndAUAU=true but cannot open PP file: " << InputPP(isRun25pp)
                       << " (mode=" << RunModeLabel(mode) << "). Skipping PP vs AuAu overlays."
                       << ANSI_RESET << "\n";
                }
          }

          if (dsPP)
          {
                analysis::RunPPvsAuAuDeliverables(*dsPP);
          }

          if (ppFallbackFile) { ppFallbackFile->Close(); delete ppFallbackFile; }
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

            const bool isCentHist = (it.path.find("_cent_") != string::npos);
            if (!ds.centSuffix.empty())
            {
              if (isCentHist && it.path.find(ds.centSuffix) == string::npos) continue;
            }

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
