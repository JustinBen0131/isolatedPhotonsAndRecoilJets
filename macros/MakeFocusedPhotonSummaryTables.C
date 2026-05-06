#include "TCanvas.h"
#include "TImage.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSystem.h"

#include <array>
#include <iostream>
#include <string>

namespace
{
const std::string kBase =
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/focusedPhotonUnfoldingQuick/pp_vs_auau020_referenceCuts";

void DrawImagePad(const std::string& path,
                  double x1,
                  double y1,
                  double x2,
                  double y2,
                  int cropTopPixels,
                  UInt_t targetW,
                  UInt_t targetH)
{
    TPad* pad = new TPad(("pad_" + path).c_str(), "", x1, y1, x2, y2);
    pad->SetFillColor(0);
    pad->SetFrameLineColor(kGray + 2);
    pad->SetFrameLineWidth(1);
    pad->SetLeftMargin(0.0);
    pad->SetRightMargin(0.0);
    pad->SetTopMargin(0.0);
    pad->SetBottomMargin(0.0);
    pad->Draw();
    pad->cd();

    TImage* img = TImage::Open(path.c_str());
    if (!img)
    {
        TLatex miss;
        miss.SetNDC();
        miss.SetTextAlign(22);
        miss.SetTextSize(0.040);
        miss.DrawLatex(0.5, 0.5, ("Missing: " + path).c_str());
        return;
    }
    if (cropTopPixels > 0)
    {
        const UInt_t width = img->GetWidth();
        const UInt_t height = img->GetHeight();
        if (height > (UInt_t)cropTopPixels + 40)
        {
            img->Crop(0, cropTopPixels, width, height - cropTopPixels);
        }
    }
    img->Scale(targetW, targetH);
    img->Draw("X");
}

void DrawOneTable(const std::string& mode,
                  const std::string& title,
                  const std::string& outPath)
{
    gROOT->SetBatch(kTRUE);

    TCanvas c(("c_summary_" + mode).c_str(), "", 2200, 1240);
    c.SetFillColor(0);
    c.SetMargin(0, 0, 0, 0);
    c.cd();

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(132);
    lat.SetTextSize(0.035);
    lat.DrawLatex(0.018, 0.972, title.c_str());

    lat.SetTextFont(62);
    lat.SetTextSize(0.021);
    lat.SetTextAlign(22);

    const std::array<std::string, 3> colLabels = {
        "Yield Overlay",
        "Response Matrix",
        "Efficiency Diagnostic"
    };
    const std::array<std::string, 2> rowLabels = {
        "pp",
        "AuAu 0-20%"
    };
    const std::array<std::string, 3> fileNames = {
        "01_photon_yield_overlay.png",
        "02_photon_response_matrix_reco_vs_truth.png",
        "03_photon_efficiency_diagnostic.png"
    };
    const std::array<std::string, 2> rowDirs = {
        "pp",
        "auau_0_20"
    };

    const double left = 0.074;
    const double right = 0.014;
    const double top = 0.115;
    const double bottom = 0.028;
    const double colGap = 0.010;
    const double rowGap = 0.045;
    const double panelW = (1.0 - left - right - 2.0 * colGap) / 3.0;
    const double panelH = (1.0 - top - bottom - rowGap) / 2.0;
    const double yBottom = bottom;
    const double yTop = bottom + panelH + rowGap;

    for (int ic = 0; ic < 3; ++ic)
    {
        const double x = left + ic * (panelW + colGap);
        lat.DrawLatex(x + 0.5 * panelW, 1.0 - top + 0.020, colLabels[ic].c_str());
    }

    lat.SetTextAngle(90);
    lat.SetTextSize(0.025);
    lat.DrawLatex(0.042, yTop + 0.5 * panelH, rowLabels[0].c_str());
    lat.DrawLatex(0.042, yBottom + 0.5 * panelH, rowLabels[1].c_str());
    lat.SetTextAngle(0);

    const std::array<int, 3> cropTop = {
        112,
        178,
        58
    };

    for (int ir = 0; ir < 2; ++ir)
    {
        const double y = (ir == 0) ? yTop : yBottom;
        for (int ic = 0; ic < 3; ++ic)
        {
            const double x = left + ic * (panelW + colGap);
            const std::string path = kBase + "/" + rowDirs[ir] + "/" + mode + "/" + fileNames[ic];
            const UInt_t targetW = (UInt_t)((panelW) * 2200.0);
            const UInt_t targetH = (UInt_t)((panelH) * 1240.0);
            DrawImagePad(path, x, y, x + panelW, y + panelH, cropTop[ic], targetW, targetH);
            c.cd();
        }
    }

    TLine sep;
    sep.SetLineColor(kGray + 1);
    sep.SetLineWidth(2);
    sep.DrawLine(left - 0.020, bottom + panelH + 0.5 * rowGap,
                 1.0 - right, bottom + panelH + 0.5 * rowGap);
    for (int ic = 1; ic < 3; ++ic)
    {
        const double x = left + ic * (panelW + colGap) - 0.5 * colGap;
        sep.DrawLine(x, bottom - 0.006, x, 1.0 - top + 0.006);
    }

    lat.SetTextFont(42);
    lat.SetTextSize(0.012);
    lat.SetTextAlign(31);
    lat.DrawLatex(0.990, 0.010, "Focused photon unfolding quick output");

    gSystem->mkdir(gSystem->DirName(outPath.c_str()), true);
    c.SaveAs(outPath.c_str());
}
} // namespace

void MakeFocusedPhotonSummaryTables()
{
    DrawOneTable("nonCorrected",
                 "Unfolded Yields & Efficiencies - no Leakage Correction",
                 kBase + "/auau_0_20/nonCorrected/00_summary_2x3_pp_over_auau_nonCorrected.png");

    DrawOneTable("leakageCorrected",
                 "Unfolded Yields & Efficiencies - Leakage Corrected",
                 kBase + "/auau_0_20/leakageCorrected/00_summary_2x3_pp_over_auau_leakageCorrected.png");

    std::cout << "[OK] " << kBase
              << "/auau_0_20/nonCorrected/00_summary_2x3_pp_over_auau_nonCorrected.png\n";
    std::cout << "[OK] " << kBase
              << "/auau_0_20/leakageCorrected/00_summary_2x3_pp_over_auau_leakageCorrected.png\n";
}
