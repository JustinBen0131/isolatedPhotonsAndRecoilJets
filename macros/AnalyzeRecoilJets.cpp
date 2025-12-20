// AnalyzeRecoilJets.cpp  —  ROOT macro to postprocess RecoilJets outputs
// Run with:  root -b -q -l AnalyzeRecoilJets.cpp
// (No function call needed; the macro calls itself at the end.)

#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TEfficiency.h>
#include <TLatex.h>

#include <iostream>
#include <iomanip>
#include <regex>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using std::cout;
using std::cerr;
using std::endl;
using std::left;
using std::right;
using std::setw;
using std::string;
using std::vector;
using std::pair;
using std::map;

// ---------------------------------------------------------------------------
// Fixed inputs/outputs (edit here only if your layout changes)
// ---------------------------------------------------------------------------
static const string kPPFile  = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp/RecoilJets_isPP_run00051282_grp008.root";
static const string kAAFile  = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/RecoilJets_isAuAu_run00068421_grp005.root";

static const string kPPPlots = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/pp/plots";
static const string kAAPlots = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auau/plots";
static const string kCombinedPlots = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/COMBINED";

// Default bins (must match what was used in production; names are parsed from hists)
static const int    kEtEdges[]   = {2,4,6,8,10,12,15,18,20,24,30};
static const size_t kNEtEdges    = sizeof(kEtEdges)/sizeof(kEtEdges[0]);

static const int    kCentEdges[] = {0,10,20,30,40,60,80};
static const size_t kNCentEdges  = sizeof(kCentEdges)/sizeof(kCentEdges[0]);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
struct SliceKey {
  int et_lo{0}, et_hi{0};
  bool hasCent{false};
  int c_lo{0}, c_hi{0};
  bool operator<(const SliceKey& o) const {
    if (et_lo != o.et_lo) return et_lo < o.et_lo;
    if (et_hi != o.et_hi) return et_hi < o.et_hi;
    if (hasCent != o.hasCent) return hasCent < o.hasCent;
    if (c_lo != o.c_lo) return c_lo < o.c_lo;
    return c_hi < o.c_hi;
  }
};

struct Counts {
  long long isoTight{0};
  long long nonIsoTight{0};
  long long isoNonTight{0};
  long long nonIsoNonTight{0};
  long long idsb_pass{0};
  long long idsb_total{0};
};

struct Eff {
  double f{0}, err{0}, lo{0}, hi{0};
  long long pass{0}, tot{0};
};

static std::regex reET(R"(.*_ET_([0-9]+)_([0-9]+))");
static std::regex reC (R"(.*_cent_([0-9]+)_([0-9]+))");

static string etTag(int lo,int hi){ return "ET_"+std::to_string(lo)+"_"+std::to_string(hi); }
static string centTag(int lo,int hi){ return "cent_"+std::to_string(lo)+"_"+std::to_string(hi); }
static string centFolder(int lo,int hi){ return "Cent_"+std::to_string(lo)+"_"+std::to_string(hi); }

static void mkdir_p(const string& path){
  if (gSystem->AccessPathName(path.c_str())) gSystem->mkdir(path.c_str(), /*recursive*/true);
}

static bool parseEtCent(const string& name, SliceKey& key){
  std::smatch m;
  // Use regex_search so we can find ET (and optional centrality) anywhere in the name
  if (!std::regex_search(name, m, reET)) return false;
  key.et_lo = std::stoi(m[1]);
  key.et_hi = std::stoi(m[2]);

  key.hasCent = false;
  if (std::regex_search(name, m, reC)) {
    key.hasCent = true;
    key.c_lo = std::stoi(m[1]);
    key.c_hi = std::stoi(m[2]);
  }
  return true;
}


static Eff eff_from_counts(long long pass, long long tot){
  Eff e; e.pass=pass; e.tot=tot;
  if (tot<=0){
    e.f=e.err=e.lo=e.hi=0;
    return e;
  }
  e.f = double(pass)/double(tot);
  // 68.27% Clopper–Pearson
  e.lo = TEfficiency::ClopperPearson((int)tot,(int)pass,0.682689,/*upper?*/false);
  e.hi = TEfficiency::ClopperPearson((int)tot,(int)pass,0.682689,/*upper?*/true);
  e.err = 0.5*(e.hi - e.lo);  // symmetric display; we still print [lo,hi]
  return e;
}

// Root directory utilities
static vector<string> triggerFolders(TFile* f){
  vector<string> out;
  TIter next(f->GetListOfKeys());
  while (auto* k = (TKey*)next()){
    if (!k->IsFolder()) continue;
    out.emplace_back(k->GetName());
  }
  std::sort(out.begin(), out.end());
  return out;
}

static void dumpAllHistEntries(const string& dataset, TFile* f){
  cout << "\n\033[1m[Inventory] Dataset: " << dataset << "\033[0m\n";
  for (const auto& trig : triggerFolders(f)){
    TDirectory* d = f->GetDirectory(trig.c_str());
    if (!d) continue;
    cout << "  Trigger: " << trig << "\n";
    TIter it(d->GetListOfKeys());
    while (auto* k = (TKey*)it()){
      TObject* obj = k->ReadObj();
      if (!obj) { cout << "    (NULL) " << k->GetName() << "\n"; continue; }
      TH1* h = dynamic_cast<TH1*>(obj);
      if (!h) { cout << "    [skip] " << k->GetName() << " (class " << obj->ClassName() << ")\n"; continue; }
      cout << "    " << std::left << setw(50) << h->GetName()
           << " entries=" << std::right << setw(10) << (Long64_t)h->GetEntries() << "\n";
    }
  }
}

// Collect sliced counters we care about
static map<SliceKey,Counts> collectCounters(TDirectory* d){
  map<SliceKey,Counts> M;
  if (!d) return M;
  TIter it(d->GetListOfKeys());
  while (auto* k = (TKey*)it()){
    TObject* obj = k->ReadObj();
    TH1* h = dynamic_cast<TH1*>(obj);
    if (!h) continue;
    const string n = h->GetName();
    SliceKey key;
    if (!parseEtCent(n, key)) continue; // skip non-sliced hists

    auto& C = M[key];
    const long long val = (Long64_t)h->GetEntries(); // 1-bin counters are written this way

    if      (n.rfind("h_isIsolated_isTight",0)==0)       C.isoTight          += val;
    else if (n.rfind("h_notIsolated_isTight",0)==0)      C.nonIsoTight       += val;
    else if (n.rfind("h_isIsolated_notTight",0)==0)      C.isoNonTight       += val;
    else if (n.rfind("h_notIsolated_notTight",0)==0)     C.nonIsoNonTight    += val;
    else if (n.rfind("h_idSB_pass",0)==0)                C.idsb_pass         += val;
    else if (n.rfind("h_idSB_total",0)==0)               C.idsb_total        += val;
  }
  return M;
}

static void saveHistAsPNG(TH1* h, const string& outPath){
  if (!h) return;
  mkdir_p(gSystem->DirName(outPath.c_str()));
  TCanvas c("c","c",1000,800);
  gStyle->SetOptStat(0);
  if      (dynamic_cast<TH2*>(h)) h->Draw("COLZ");
  else if (dynamic_cast<TH3*>(h)) h->Draw("LEGO2Z");
  else                            h->Draw();
  c.SaveAs(outPath.c_str());
}

static void dumpFolderPNGs(TDirectory* d, const string& outDir){
  mkdir_p(outDir);
  TIter it(d->GetListOfKeys());
  while (auto* k = (TKey*)it()){
    TObject* obj = k->ReadObj();
    TH1* h = dynamic_cast<TH1*>(obj);
    if (!h) continue;
    const string out = outDir + "/" + string(h->GetName()) + ".png";
    saveHistAsPNG(h, out);
  }
}

// Filter and dump per-slice histograms into requested ET/Cent trees
static void dumpSlicePNGs(TDirectory* d, const string& baseOut, bool isAuAu){
  TIter it(d->GetListOfKeys());
  while (auto* k = (TKey*)it()){
    TObject* obj = k->ReadObj();
    TH1* h = dynamic_cast<TH1*>(obj);
    if (!h) continue;
    const string name = h->GetName();
    SliceKey key;
    if (!parseEtCent(name, key)) continue;

    if (!isAuAu){
      // pp: /<TRIG>/ET_lo_hi/
      const string dir = baseOut + "/" + etTag(key.et_lo,key.et_hi);
      saveHistAsPNG(h, dir + "/" + name + ".png");
    } else {
      if (!key.hasCent) continue; // AuAu slice images go under Cent_x_y/ET_lo_hi
      const string dir = baseOut + "/" + centFolder(key.c_lo,key.c_hi) + "/" + etTag(key.et_lo,key.et_hi);
      saveHistAsPNG(h, dir + "/" + name + ".png");
    }
  }
}

// Make a simple TGraphErrors from (x,y,err) triplets
static TGraphErrors* makeGraph(const vector<double>& x, const vector<double>& y, const vector<double>& ey, int color, int marker){
  auto* g = new TGraphErrors((int)x.size());
  for (size_t i=0;i<x.size();++i){
    g->SetPoint((int)i, x[i], y[i]);
    g->SetPointError((int)i, 0.0, (i<ey.size()?ey[i]:0.0));
  }
  g->SetLineColor(color); g->SetMarkerColor(color);
  g->SetMarkerStyle(marker); g->SetMarkerSize(1.2); g->SetLineWidth(2);
  return g;
}

// Draw one graph with axes and optional overlays, save png
static void drawAndSaveGraphs(const string& title, const string& ytitle,
                              const string& outPng, const vector<TGraphErrors*>& graphs,
                              const vector<string>& labels,
                              double ymin=0.0, double ymax=1.0)
{
  mkdir_p(gSystem->DirName(outPng.c_str()));
  TCanvas c("c","c",1100,850);
  gStyle->SetOptStat(0);

  // Determine x range from ET bins
  double xmin = kEtEdges[0], xmax = kEtEdges[kNEtEdges-1];

  // Draw a dummy frame
  TH1F frame("frame","",100,xmin-0.5,xmax+0.5);
  frame.GetXaxis()->SetTitle("E_{T}^{#gamma} [GeV]");
  frame.GetYaxis()->SetTitle(ytitle.c_str());
  frame.SetMinimum(ymin); frame.SetMaximum(ymax);
  frame.Draw();

  // Title annotation
  TLatex lat; lat.SetNDC(); lat.SetTextSize(0.035);
  lat.DrawLatex(0.14,0.93,Form("#bf{%s}", title.c_str()));

  // Legend
  TLegend leg(0.62,0.70,0.92,0.92);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);

  for (size_t i=0;i<graphs.size();++i){
    if (!graphs[i]) continue;
    graphs[i]->Draw("P SAME");
    if (i<labels.size()) leg.AddEntry(graphs[i], labels[i].c_str(), "pe");
  }
  leg.Draw();
  c.SaveAs(outPng.c_str());
}

// ---------------------------------------------------------------------------
// Pretty terminal tables
// ---------------------------------------------------------------------------
static void printPPTable(const string& trig, const map<SliceKey,Counts>& M){
  cout << "\n\033[1m[pp] Trigger: " << trig << " — per‑ET counts and ID‑SB fraction\033[0m\n";
  cout << " ET bin   |  iso^tight | !iso^tight | iso^!tight | !iso^!tight |  SB pass/total  | f_IDSB  ±  σ   [lo,hi]\n";
  cout << "----------+------------+------------+------------+-------------+-----------------+-------------------------\n";
  for (size_t i=0;i+1<kNEtEdges;++i){
    SliceKey key; key.et_lo=kEtEdges[i]; key.et_hi=kEtEdges[i+1]; key.hasCent=false;
    auto it = M.find(key);
    Counts C;
    if (it!=M.end()) C = it->second;
    Eff e = eff_from_counts(C.idsb_pass, C.idsb_total);
    printf(" %2d–%-3d  | %10lld | %10lld | %10lld | %11lld | %7lld/%-7lld | %6.3f ± %5.3f  [%5.3f,%5.3f]\n",
           key.et_lo, key.et_hi,
           C.isoTight, C.nonIsoTight, C.isoNonTight, C.nonIsoNonTight,
           C.idsb_pass, C.idsb_total, e.f, e.err, e.lo, e.hi);
  }
}

static void printAATables(const string& trig, const map<SliceKey,Counts>& M){
  cout << "\n\033[1m[Au+Au] Trigger: " << trig << " — per‑centrality, per‑ET counts and ID‑SB fraction\033[0m\n";
  for (size_t jc=0;jc+1<kNCentEdges;++jc){
    const int clo=kCentEdges[jc], chi=kCentEdges[jc+1];
    cout << "\n  \033[1mCentrality: " << clo << "–" << chi << " %\033[0m\n";
    cout << "  ET bin  |  iso^tight | !iso^tight | iso^!tight | !iso^!tight |  SB pass/total  | f_IDSB  ±  σ   [lo,hi]\n";
    cout << "  --------+------------+------------+------------+-------------+-----------------+-------------------------\n";
    for (size_t i=0;i+1<kNEtEdges;++i){
      SliceKey key; key.et_lo=kEtEdges[i]; key.et_hi=kEtEdges[i+1]; key.hasCent=true; key.c_lo=clo; key.c_hi=chi;
      auto it = M.find(key);
      Counts C;
      if (it!=M.end()) C = it->second;
      Eff e = eff_from_counts(C.idsb_pass, C.idsb_total);
      printf("  %2d–%-3d | %10lld | %10lld | %10lld | %11lld | %7lld/%-7lld | %6.3f ± %5.3f  [%5.3f,%5.3f]\n",
            key.et_lo, key.et_hi,
            C.isoTight, C.nonIsoTight, C.isoNonTight, C.nonIsoNonTight,
            C.idsb_pass, C.idsb_total, e.f, e.err, e.lo, e.hi);
    }
  }
}

static void printRCTable(const string& caption,
                         const map<SliceKey,Counts>& M_AA,
                         const map<SliceKey,Counts>& M_PP)
{
  cout << "\n\033[1m[r_c] " << caption << "\033[0m\n";
  cout << "  Cent     |  ET bin  |   f_IDSB^pp ± σ     |  f_IDSB^AuAu ± σ    |   r_c ± δr\n";
  cout << "  ---------+----------+----------------------+----------------------+----------------\n";
  for (size_t jc=0;jc+1<kNCentEdges;++jc){
    const int clo=kCentEdges[jc], chi=kCentEdges[jc+1];
    for (size_t i=0;i+1<kNEtEdges;++i){
      SliceKey ka; ka.et_lo=kEtEdges[i]; ka.et_hi=kEtEdges[i+1]; ka.hasCent=true; ka.c_lo=clo; ka.c_hi=chi;
      SliceKey kp; kp.et_lo=kEtEdges[i]; kp.et_hi=kEtEdges[i+1]; kp.hasCent=false;

      Counts Caa,Cpp;
      auto ita=M_AA.find(ka); if (ita!=M_AA.end()) Caa = ita->second;
      auto itp=M_PP.find(kp); if (itp!=M_PP.end()) Cpp = itp->second;

      Eff eAA = eff_from_counts(Caa.idsb_pass, Caa.idsb_total);
      Eff ePP = eff_from_counts(Cpp.idsb_pass, Cpp.idsb_total);

      double rc=0, dr=0;
      if (ePP.tot>0 && ePP.f>0 && eAA.tot>0){
        rc = eAA.f / ePP.f;
        // error propagation (Wald)
        const double sPP = (ePP.tot>0? std::sqrt(ePP.f*(1-ePP.f)/ePP.tot) : 0);
        const double sAA = (eAA.tot>0? std::sqrt(eAA.f*(1-eAA.f)/eAA.tot) : 0);
        if (eAA.f>0 && ePP.f>0) dr = rc*std::sqrt( (sAA/eAA.f)*(sAA/eAA.f) + (sPP/ePP.f)*(sPP/ePP.f) );
      }

      printf("  %2d–%-3d  | %2d–%-3d  |  %6.3f ± %5.3f     |  %6.3f ± %5.3f     |  %6.3f ± %5.3f\n",
             clo,chi, kEtEdges[i],kEtEdges[i+1],
             ePP.f, (ePP.tot? std::sqrt(ePP.f*(1-ePP.f)/ePP.tot):0.0),
             eAA.f, (eAA.tot? std::sqrt(eAA.f*(1-eAA.f)/eAA.tot):0.0),
             rc,  dr);
    }
  }
}

// ---------------------------------------------------------------------------
// Plot builders (f_IDSB and r_c)
// ---------------------------------------------------------------------------
static void plot_fIDSB_pp(const string& trig,
                          const map<SliceKey,Counts>& M,
                          const string& outDir)
{
  vector<double> x,y,ey;
  for (size_t i=0;i+1<kNEtEdges;++i){
    SliceKey k; k.et_lo=kEtEdges[i]; k.et_hi=kEtEdges[i+1]; k.hasCent=false;
    auto it=M.find(k); if (it==M.end()) continue;
    Eff e = eff_from_counts(it->second.idsb_pass, it->second.idsb_total);
    if (e.tot<=0) continue;
    x.push_back(0.5*(k.et_lo+k.et_hi)); y.push_back(e.f); ey.push_back(e.err);
  }
  auto* g = makeGraph(x,y,ey,kBlue+1,20);
  const string png = outDir + "/f_IDSB_vs_ET.png";
  drawAndSaveGraphs(trig+"  (pp)", "f_{ID-SB}", png, {g}, {trig}, 0.0, 1.05);
}

static void plot_fIDSB_aa(const string& trig,
                          const map<SliceKey,Counts>& M,
                          const string& trigOut)
{
  // Per‑centrality graphs + overlay
  vector<TGraphErrors*> overlay;
  vector<string> labels;
  int colors[] = {kRed+1,kBlue+1,kGreen+3,kMagenta+1,kOrange+7,kTeal+2};
  int markers[] = {20,21,22,23,33,34};

  for (size_t jc=0;jc+1<kNCentEdges;++jc){
    const int clo=kCentEdges[jc], chi=kCentEdges[jc+1];
    vector<double> x,y,ey;
    for (size_t i=0;i+1<kNEtEdges;++i){
      SliceKey k; k.et_lo=kEtEdges[i]; k.et_hi=kEtEdges[i+1]; k.hasCent=true; k.c_lo=clo; k.c_hi=chi;
      auto it=M.find(k); if (it==M.end()) continue;
      Eff e = eff_from_counts(it->second.idsb_pass, it->second.idsb_total);
      if (e.tot<=0) continue;
      x.push_back(0.5*(k.et_lo+k.et_hi)); y.push_back(e.f); ey.push_back(e.err);
    }
    if (x.empty()) continue;

    const string centDir = trigOut + "/" + centFolder(clo,chi);
    mkdir_p(centDir);
    auto* g = makeGraph(x,y,ey,colors[jc%6],markers[jc%6]);
    drawAndSaveGraphs(Form("%s  (Au+Au, %d–%d%%)", trig.c_str(), clo,chi),
                      "f_{ID-SB}", centDir+"/f_IDSB_vs_ET.png", {g}, {Form("%d–%d%%",clo,chi)}, 0.0, 1.05);
    overlay.push_back(g);
    labels.emplace_back(Form("%d–%d%%",clo,chi));
  }
  if (!overlay.empty()){
    drawAndSaveGraphs(trig+"  (Au+Au)", "f_{ID-SB}",
                      trigOut+"/f_IDSB_vs_ET_byCentrality.png", overlay, labels, 0.0, 1.05);
  }
}

static void plot_rc_overlay(const map<SliceKey,Counts>& M_AA,
                            const map<SliceKey,Counts>& M_PP,
                            const string& outCombined,
                            const string& overlayTitle)
{
  mkdir_p(outCombined);
  vector<TGraphErrors*> over; vector<string> labels;
  int colors[]  ={kRed+1,kBlue+1,kGreen+3,kMagenta+1,kOrange+7,kTeal+2};
  int markers[] ={20,21,22,23,33,34};

  for (size_t jc=0;jc+1<kNCentEdges;++jc){
    const int clo=kCentEdges[jc], chi=kCentEdges[jc+1];
    vector<double> x,y,ey;
    for (size_t i=0;i+1<kNEtEdges;++i){
      SliceKey ka; ka.et_lo=kEtEdges[i]; ka.et_hi=kEtEdges[i+1]; ka.hasCent=true; ka.c_lo=clo; ka.c_hi=chi;
      SliceKey kp; kp.et_lo=kEtEdges[i]; kp.et_hi=kEtEdges[i+1]; kp.hasCent=false;
      auto ita=M_AA.find(ka); auto itp=M_PP.find(kp);
      if (ita==M_AA.end() || itp==M_PP.end()) continue;

      Eff eAA = eff_from_counts(ita->second.idsb_pass, ita->second.idsb_total);
      Eff ePP = eff_from_counts(itp->second.idsb_pass, itp->second.idsb_total);
      if (ePP.tot<=0 || eAA.tot<=0 || ePP.f<=0) continue;

      const double rc = eAA.f / ePP.f;
      const double sPP = std::sqrt(ePP.f*(1-ePP.f)/ePP.tot);
      const double sAA = std::sqrt(eAA.f*(1-eAA.f)/eAA.tot);
      const double dr = rc*std::sqrt( (sAA/eAA.f)*(sAA/eAA.f) + (sPP/ePP.f)*(sPP/ePP.f) );

      x.push_back(0.5*(ka.et_lo+ka.et_hi)); y.push_back(rc); ey.push_back(dr);
    }
    if (x.empty()) continue;
    auto* g = makeGraph(x,y,ey,colors[jc%6],markers[jc%6]);
    over.push_back(g); labels.emplace_back(Form("%d–%d%%",clo,chi));

    // Individual centrality canvas too
    drawAndSaveGraphs(Form("r_{c}(c,ET)  (Cent %d–%d%%)",clo,chi),
                      "r_{c} = f^{AuAu}_{ID-SB} / f^{pp}_{ID-SB}",
                      outCombined + Form("/rc_vs_ET_Cent_%d_%d.png",clo,chi),
                      {g},{Form("%d–%d%%",clo,chi)}, 0.0, 2.0);
  }

  if (!over.empty()){
    drawAndSaveGraphs(overlayTitle,
                      "r_{c} = f^{AuAu}_{ID-SB} / f^{pp}_{ID-SB}",
                      outCombined + "/rc_vs_ET_overlay.png",
                      over, labels, 0.0, 2.0);
  }
}

// ---------------------------------------------------------------------------
// Main driver
// ---------------------------------------------------------------------------
static void AnalyzeRecoilJets()
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  // -------- Open files (fail loudly if missing) --------
  std::unique_ptr<TFile> fPP( TFile::Open(kPPFile.c_str(),"READ") );
  if (!fPP || fPP->IsZombie()){ cerr << "[ERROR] Cannot open " << kPPFile << "\n"; return; }

  std::unique_ptr<TFile> fAA( TFile::Open(kAAFile.c_str(),"READ") );
  if (!fAA || fAA->IsZombie()){ cerr << "[ERROR] Cannot open " << kAAFile << "\n"; return; }

  // -------- Full inventory with entries --------
  dumpAllHistEntries("pp",   fPP.get());
  dumpAllHistEntries("Au+Au",fAA.get());

  // -------- Process each dataset trigger-by-trigger --------
  // (We assume exactly one relevant trigger per file, but handle >1 robustly.)
  auto trigPP = triggerFolders(fPP.get());
  auto trigAA = triggerFolders(fAA.get());

  if (trigPP.empty()) cerr << "[WARN] No top-level trigger folders found in pp file.\n";
  if (trigAA.empty()) cerr << "[WARN] No top-level trigger folders found in Au+Au file.\n";

  // Keep a single (last seen) map to use for the combined rc plots;
  // if multiple triggers exist, we use the first for each dataset.
  map<SliceKey,Counts> firstPP, firstAA;
  string firstTrigPP, firstTrigAA;

  // ---------- pp loop ----------
  for (size_t it=0; it<trigPP.size(); ++it){
    const string trig = trigPP[it];
    TDirectory* d = fPP->GetDirectory(trig.c_str());
    if (!d){ cerr << "[WARN] Missing trigger folder in pp: " << trig << "\n"; continue; }

    // Counters we need
    auto M = collectCounters(d);

    // Print per-ET table
    printPPTable(trig, M);

    // Plots: dump all hists and per-ET folders
    const string trigOut = kPPPlots + "/" + trig;
    mkdir_p(trigOut);
    dumpFolderPNGs(d, trigOut + "/allHists");
    dumpSlicePNGs(d, trigOut, /*isAuAu=*/false);

    // f_IDSB vs ET
    plot_fIDSB_pp(trig, M, trigOut);

    if (firstPP.empty()){ firstPP = M; firstTrigPP = trig; }
  }

  // ---------- Au+Au loop ----------
  for (size_t it=0; it<trigAA.size(); ++it){
    const string trig = trigAA[it];
    TDirectory* d = fAA->GetDirectory(trig.c_str());
    if (!d){ cerr << "[WARN] Missing trigger folder in Au+Au: " << trig << "\n"; continue; }

    auto M = collectCounters(d);

    // Print per-cent, per-ET tables
    printAATables(trig, M);

    // Plots
    const string trigOut = kAAPlots + "/" + trig;
    mkdir_p(trigOut);
    dumpFolderPNGs(d, trigOut + "/allHists");
    dumpSlicePNGs(d, trigOut, /*isAuAu=*/true);

    // f_IDSB per-centrality + overlay
    plot_fIDSB_aa(trig, M, trigOut);

    if (firstAA.empty()){ firstAA = M; firstTrigAA = trig; }
  }

  // ---------- r_c tables + plots (combined) ----------
  if (!firstAA.empty() && !firstPP.empty()){
    printRCTable("Using first trigger in each file "
                 "(" + (firstTrigAA.empty()?"AA?":firstTrigAA) + " vs "
                      + (firstTrigPP.empty()?"PP?":firstTrigPP) + ")",
                 firstAA, firstPP);

    plot_rc_overlay(firstAA, firstPP, kCombinedPlots,
                    "r_{c}(c,ET)  —  " + (firstTrigAA.empty()?"AA?":firstTrigAA)
                    + "  vs  " + (firstTrigPP.empty()?"PP?":firstTrigPP));
  } else {
    cerr << "[WARN] Skipping r_c plots (missing counters in pp or Au+Au).\n";
  }

  cout << "\n\033[1;32mAll done.\033[0m\n";
}

// ---------------------------------------------------------------------------
// Self-invoke when run as a ROOT macro
// ---------------------------------------------------------------------------
void __AnalyzeRecoilJets_autorun__(){ AnalyzeRecoilJets(); }

