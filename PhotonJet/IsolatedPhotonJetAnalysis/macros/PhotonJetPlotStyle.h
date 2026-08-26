#pragma once

// Canonical ROOT-side presentation of a generated PhotonJetPlotContractV1.
// Physics labels are supplied only by the generated annotation header; this
// helper intentionally has no string arguments for cuts, samples, or status.

#include <TLatex.h>
#include <TPad.h>
#include <TStyle.h>

#include <cstddef>
#include <stdexcept>

namespace PhotonJetPlotStyle
{
  inline void Apply()
  {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetPadColor(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
  }

  template <class GeneratedAnnotation>
  inline void DrawGeneratedHeader(const GeneratedAnnotation& annotation)
  {
    if (!gPad)
    {
      throw std::runtime_error("DrawGeneratedHeader requires an active ROOT pad");
    }
    // Individual plots keep annotations inside the frame. A large external
    // top margin indicates the obsolete whitespace-header layout.
    if (gPad->GetTopMargin() > 0.12)
    {
      throw std::runtime_error("canonical individual plots forbid an external header band");
    }

    TLatex text;
    text.SetNDC(true);
    text.SetTextFont(42);
    text.SetTextAlign(13);
    text.SetTextSize(0.044);
    text.DrawLatex(0.15, 0.88, annotation.kExperimentLabel.data());

    text.SetTextSize(0.030);
    for (std::size_t index = 0; index < annotation.kLines.size(); ++index)
    {
      const double x = 0.15;
      const double y = 0.825 - 0.045 * index;
      text.DrawLatex(x, y, annotation.kLines[index].data());
    }
  }
}  // namespace PhotonJetPlotStyle
