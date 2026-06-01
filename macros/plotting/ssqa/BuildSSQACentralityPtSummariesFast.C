#include <string>
#include <vector>

#include "TCanvas.h"
#include "TImage.h"
#include "TLatex.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TSystem.h"

namespace {

std::string JoinPathFast(const std::string& a, const std::string& b)
{
    if (a.empty()) return b;
    if (b.empty()) return a;
    if (a.back() == '/') return a + b;
    return a + "/" + b;
}

void EnsureDirFast(const std::string& dir)
{
    if (!dir.empty()) gSystem->mkdir(dir.c_str(), true);
}

void DrawMissingPadFast(const std::string& titleLine)
{
    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextAlign(22);
    t.SetTextSize(0.080);
    t.DrawLatex(0.50, 0.55, "MISSING");

    t.SetTextSize(0.050);
    t.DrawLatex(0.50, 0.42, titleLine.c_str());
}

} // namespace

int BuildSSQACentralityPtSummariesFast(const char* variantFolder = "baseVariant")
{
    const std::string ssQADir =
        "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/"
        "jetMinPt7_7pi_8_vz30_isoR30_isSliding/"
        "photon_10_plus_MBD_NS_geq_2_vtx_lt_150/SS_QA";

    const std::string centPtSummaryDir = JoinPathFast(ssQADir, "centralityXpTsummaries");
    const std::string variantSummaryDir = JoinPathFast(centPtSummaryDir, variantFolder);
    EnsureDirFast(variantSummaryDir);

    const std::vector<std::string> tags = {"inclusive", "pre", "tight", "nonTight"};
    const std::vector<std::string> vars = {
        "weta", "wphi", "weta35", "wphi53", "e11e33", "et1", "e32e35"
    };
    const std::vector<std::string> summaryCentFolders = {"0_20", "20_50", "50_80"};
    const std::vector<std::string> summaryCentLabels  = {"0-20%", "20-50%", "50-80%"};
    const std::vector<std::string> summaryPtFolders   = {"pT_10_12", "pT_14_16", "pT_18_20", "pT_20_35"};
    const std::vector<std::string> summaryPtLabels    = {"10-12 GeV", "14-16 GeV", "18-20 GeV", "20-35 GeV"};

    for (const auto& tag : tags)
    {
        const std::string tagSummaryDir = JoinPathFast(variantSummaryDir, tag);
        EnsureDirFast(tagSummaryDir);

        for (const auto& var : vars)
        {
            const std::string cSummaryName = "c_ssQA_centXpT_" + std::string(variantFolder) + "_" + tag + "_" + var + "_fast";
            TCanvas cSummary(cSummaryName.c_str(), cSummaryName.c_str(), 5200, 3000);
            cSummary.cd();
            cSummary.SetFillColor(0);
            cSummary.SetFrameFillColor(0);

            const double leftSummary   = 0.055;
            const double rightSummary  = 0.010;
            const double topSummary    = 0.080;
            const double bottomSummary = 0.018;
            const double xGapSummary   = 0.004;
            const double yGapSummary   = 0.010;

            const int nSummaryRows = (int)summaryCentFolders.size();
            const int nSummaryCols = (int)summaryPtFolders.size();
            const double rawCellW = (1.0 - leftSummary - rightSummary) / nSummaryCols;
            const double rawCellH = (1.0 - topSummary - bottomSummary) / nSummaryRows;

            std::vector<TImage*> keepSummaryImages;
            bool anySummaryCell = false;

            for (int irow = 0; irow < nSummaryRows; ++irow)
            {
                for (int icol = 0; icol < nSummaryCols; ++icol)
                {
                    const double x1 = leftSummary + icol * rawCellW + 0.5 * xGapSummary;
                    const double x2 = leftSummary + (icol + 1) * rawCellW - 0.5 * xGapSummary;
                    const double y2 = 1.0 - topSummary - irow * rawCellH - 0.5 * yGapSummary;
                    const double y1 = 1.0 - topSummary - (irow + 1) * rawCellH + 0.5 * yGapSummary;

                    const std::string padName =
                        cSummaryName + "_pad_r" + std::to_string(irow) + "_c" + std::to_string(icol);

                    TPad* cellPad = new TPad(padName.c_str(), padName.c_str(), x1, y1, x2, y2);
                    cellPad->SetLeftMargin(0.0);
                    cellPad->SetRightMargin(0.0);
                    cellPad->SetBottomMargin(0.0);
                    cellPad->SetTopMargin(0.0);
                    cellPad->SetBorderMode(0);
                    cellPad->SetBorderSize(0);
                    cellPad->SetFillColor(0);
                    cellPad->SetFrameFillColor(0);
                    cellPad->Draw();
                    cellPad->cd();

                    const std::string cellPng = JoinPathFast(
                        JoinPathFast(
                            JoinPathFast(
                                JoinPathFast(
                                    JoinPathFast(
                                        JoinPathFast(ssQADir, summaryCentFolders[irow]),
                                        summaryPtFolders[icol]),
                                    "inclusiveBACKGROUNDphotonJetSIGNALoverlay"),
                                variantFolder),
                            "perVariable"),
                        JoinPathFast(tag, var + ".png"));

                    if (!gSystem->AccessPathName(cellPng.c_str()))
                    {
                        TImage* img = TImage::Open(cellPng.c_str());
                        if (img)
                        {
                            img->Draw();
                            TPaveText mask(0.60, 0.72, 0.998, 0.998, "NDC");
                            mask.SetBorderSize(0);
                            mask.SetFillColor(kWhite);
                            mask.SetFillStyle(1001);
                            mask.SetLineColor(kWhite);
                            mask.Draw();
                            keepSummaryImages.push_back(img);
                            anySummaryCell = true;
                        }
                        else
                        {
                            DrawMissingPadFast(summaryCentLabels[irow] + ", " + summaryPtLabels[icol]);
                        }
                    }
                    else
                    {
                        DrawMissingPadFast(summaryCentLabels[irow] + ", " + summaryPtLabels[icol]);
                    }

                    cSummary.cd();
                }
            }

            TLatex tCol;
            tCol.SetNDC(true);
            tCol.SetTextFont(42);
            tCol.SetTextAlign(22);
            tCol.SetTextSize(0.032);
            for (int icol = 0; icol < nSummaryCols; ++icol)
            {
                const double x = leftSummary + (icol + 0.5) * rawCellW;
                const double y = 1.0 - 0.58 * topSummary;
                tCol.DrawLatex(x, y, summaryPtLabels[icol].c_str());
            }

            TLatex tRow;
            tRow.SetNDC(true);
            tRow.SetTextFont(42);
            tRow.SetTextAlign(22);
            tRow.SetTextAngle(270);
            tRow.SetTextSize(0.033);
            for (int irow = 0; irow < nSummaryRows; ++irow)
            {
                const double x = 0.50 * leftSummary;
                const double y = 1.0 - topSummary - (irow + 0.5) * rawCellH;
                tRow.DrawLatex(x, y, summaryCentLabels[irow].c_str());
            }

            if (anySummaryCell)
            {
                cSummary.SaveAs(JoinPathFast(tagSummaryDir, var + ".png").c_str());
            }

            for (TImage* img : keepSummaryImages) delete img;
        }
    }

    return 0;
}
