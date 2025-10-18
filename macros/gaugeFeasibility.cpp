// gaugeFeasibilityGrid.C
// 3×3 centrality grid: Gaussian UE tails (PPG04 randomized ηφ, iterative)
// One vertical iso line at E_cut; shades left/right; prints P(E>cut) per bin.

#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <string>

static constexpr double kCut   = 4.0;          // GeV, draw line at E_T^iso = 4
static constexpr double kXmin  = -10.0;        // GeV
static constexpr double kXmax  =  30.0;        // GeV
static constexpr double kAlpha =  0.30;        // shading alpha

// Normalized Gaussian pdf
double pdf(double *x, double *p){
  const double mu=p[0], s=p[1], z=(x[0]-mu)/s;
  const double norm = 1.0 / ( TMath::Sqrt(2.0*TMath::Pi()) * s );
  return norm * TMath::Exp(-0.5*z*z);
}

// Shade area under f from xa to xb
void Shade(const TF1* f, double xa, double xb, int col, float a=0.3, int N=300){
  const double dx=(xb-xa)/N;
  for(int i=0;i<N;++i){
    const double x1=xa+i*dx, x2=x1+dx;
    const double y = TMath::Max(f->Eval(x1), f->Eval(x2));
      auto b = new TBox(x1,0,x2,y);
      b->SetFillColorAlpha(col,a);
      b->SetLineColorAlpha(col,0);
      b->Draw("f");

  }
}

void gaugeFeasibility(){
  gStyle->SetOptStat(0);
    
    // Use the function/pad title and make it prominent
    gStyle->SetOptTitle(1);         // ensure titles are drawn
    gStyle->SetTitleFont(42, "t");  // title font
    gStyle->SetTitleFontSize(0.06); // title size (try 0.06–0.08)
    gStyle->SetTitleX(0.5);         // centered
    gStyle->SetTitleY(0.97);        // near the top edge
    gStyle->SetTitleAlign(23);      // center horizontally & vertically
    gStyle->SetTitleBorderSize(0);  // no border
    gStyle->SetTitleFillColor(0);   // transparent background


  // --- centrality labels (9 bins => 3×3 grid) ---
  std::vector<std::string> lab = {
    "0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"
  };

  // --- μ,σ from PPG04 Table 4 (Randomized ηφ, Iterative) ---
  std::vector<double> mu  = { 0.31,  /*5-10*/ 0.19, /*10-20*/ 0.09, 0.05,0.016,-0.0022,-5.5e-5,0.0019,0.0067 };
  std::vector<double> sig = { 4.80,  /*5-10*/ 4.40, /*10-20*/ 3.80, 3.10,2.60,2.10,1.70,1.40,1.20 };

  // --- sanity: if any σ<=0, skip drawing tails but still draw frame ---
  auto Pright = [](double cut,double mu,double s){
    return (s>0) ? 0.5*TMath::Erfc( (cut-mu)/(TMath::Sqrt2()*s) ) : 0.0;
  };

  // --- canvas 3×3 ---
  TCanvas c("c","UE tail gauge (PPG04 randomized η#varphi, iterative)", 1200, 1000);
  c.Divide(3,3,0.002,0.002);

    std::cout.setf(std::ios::scientific);
    std::cout.precision(2);
    std::cout << "\nCentrality  mu[GeV]  sigma[GeV]  P(E>cut)  P(E<cut)  Right/Left\n";


  for (int k=0;k<9;++k){
    c.cd(k+1);

      auto f = new TF1(Form("f_%d",k), pdf, kXmin, kXmax, 2);
      f->SetParameters(mu[k], sig[k]);
      f->SetLineWidth(3);
      f->SetTitle(Form("%s;UE cone E_{T}^{iso}  [GeV];Probability density", lab[k].c_str()));
      f->Draw();

      auto L = new TLine(kCut, 0, kCut, 1.2*f->Eval(kCut));
      L->SetLineColor(kBlack); L->SetLineWidth(4); L->SetLineStyle(7); L->Draw();

    double pR = Pright(kCut, mu[k], sig[k]);
    double pL = 1.0 - pR;
    double r  = (pL>0 ? pR/pL : 0.0);

    if (sig[k]>0){
        Shade(f, kXmin, kCut, kAzure+1, kAlpha);    // left (good)
        Shade(f, kCut,  kXmax, kGreen+2, kAlpha);   // right (UE > cut)
    }


      // Top-right labels in NDC (independent of axis ranges)
      auto tx2 = new TLatex();
      tx2->SetNDC(true);                 // use normalized device coords
      tx2->SetTextAlign(33);             // right-top alignment
      tx2->SetTextFont(42);
      tx2->SetTextSize(0.052);
      tx2->DrawLatex(0.83, 0.87,
                     Form("#mu=% .3g, #sigma=% .3g   E_{cut}=%.1f", mu[k], sig[k], kCut));


      auto tx3a = new TLatex(); tx3a->SetNDC(true); tx3a->SetTextAlign(33);
      tx3a->SetTextFont(42);    tx3a->SetTextSize(0.05);
      tx3a->DrawLatex(0.865, 0.60, Form("P(E>cut) = % .4g", pR));

      auto tx3b = new TLatex(); tx3b->SetNDC(true); tx3b->SetTextAlign(33);
      tx3b->SetTextFont(42);    tx3b->SetTextSize(0.05);
      tx3b->DrawLatex(0.86, 0.52, Form("P(E<cut) = % .4g", pL));

      auto tx3c = new TLatex(); tx3c->SetNDC(true); tx3c->SetTextAlign(33);
      tx3c->SetTextFont(42);    tx3c->SetTextSize(0.05);
      tx3c->DrawLatex(0.8, 0.44, Form("Ratio = % .3g", r));


      std::cout << lab[k] << "   "
                << mu[k] << "      " << sig[k] << "        "
                << pR << "     " << pL << "     " << r << "\n";

      gPad->Modified();
      gPad->Update();
  }

    c.Modified();      // force canvas to register all drawings
    c.Update();        // flush to the off-screen buffer
    c.SaveAs("rc_grid.png");
    std::cout << "\nSaved rc_grid.png\n\n";
}

