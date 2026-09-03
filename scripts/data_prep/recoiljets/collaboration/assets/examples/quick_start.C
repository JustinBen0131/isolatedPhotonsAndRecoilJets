#include <TCanvas.h>
#include <TChain.h>
#include <TError.h>
#include <TH1D.h>
#include <TMath.h>

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <tuple>

struct QuickPhotonLeader {
  double photon_et = -1;
  int encounter_ordinal = -1;
  ULong64_t candidate_id_hi = 0;
  ULong64_t candidate_id_lo = 0;
};

static bool quick_start_add_files(TChain &chain, const char *path)
{
  const std::string value = path ? path : "";
  if (value.size() > 5 && value.substr(value.size() - 5) == ".root")
    return chain.Add(value.c_str()) > 0;
  std::ifstream input(value);
  std::string filename;
  while (std::getline(input, filename))
    if (!filename.empty() && filename[0] != '#') chain.Add(filename.c_str());
  return input.good() || input.eof();
}

void quick_start(const char *product = "../data/auau/files.txt",
                 const char *system = "auau")
{
  const std::string collision_system = system ? system : "";
  if (collision_system != "pp" && collision_system != "auau") {
    Error("quick_start", "system must be pp or auau"); return;
  }

  TChain photons("photons");
  TChain pairs("photonJets");
  if (!quick_start_add_files(photons, product) ||
      !quick_start_add_files(pairs, product)) {
    Error("quick_start", "Cannot load photons and photonJets from %s", product); return;
  }

  Int_t source_file_index = -1;
  ULong64_t event_id_hi = 0, event_id_lo = 0;
  ULong64_t candidate_id_hi = 0, candidate_id_lo = 0;
  Int_t encounter_ordinal = -1, tight = 0, isolated_r04 = 0;
  double photon_et = 0, photon_eta = 0, vertex_z = 0, centrality = 0;
  photons.SetBranchAddress("source_file_index", &source_file_index);
  photons.SetBranchAddress("event_id_hi", &event_id_hi);
  photons.SetBranchAddress("event_id_lo", &event_id_lo);
  photons.SetBranchAddress("candidate_id_hi", &candidate_id_hi);
  photons.SetBranchAddress("candidate_id_lo", &candidate_id_lo);
  photons.SetBranchAddress("photon_encounter_ordinal", &encounter_ordinal);
  photons.SetBranchAddress("photon_et", &photon_et);
  photons.SetBranchAddress("photon_eta", &photon_eta);
  photons.SetBranchAddress("vertex_z", &vertex_z);
  photons.SetBranchAddress("centrality", &centrality);
  photons.SetBranchAddress("bdt_is_tight", &tight);
  photons.SetBranchAddress("iso_r04_pass", &isolated_r04);

  using EventKey = std::tuple<Int_t, Int_t, ULong64_t, ULong64_t>;
  std::map<EventKey, QuickPhotonLeader> leaders;
  for (Long64_t i = 0; i < photons.GetEntries(); ++i) {
    photons.GetEntry(i);
    if (photon_et < 15 || photon_et >= 35 || std::abs(photon_eta) >= 0.7) continue;
    if (collision_system == "pp" && std::abs(vertex_z) >= 30) continue;
    if (collision_system == "auau" &&
        (std::abs(vertex_z) >= 10 || centrality < 0 || centrality >= 20)) continue;
    if (!tight || !isolated_r04) continue;

    QuickPhotonLeader &leader = leaders[
      {photons.GetTreeNumber(), source_file_index, event_id_hi, event_id_lo}
    ];
    const bool precedes = leader.photon_et < 0 || photon_et > leader.photon_et ||
      (photon_et == leader.photon_et && encounter_ordinal < leader.encounter_ordinal);
    if (precedes)
      leader = {photon_et, encounter_ordinal, candidate_id_hi, candidate_id_lo};
  }

  double jet_radius = 0, jet_pt = 0, jet_eta = 0;
  double delta_phi = 0, xjgamma = 0;
  pairs.SetBranchAddress("source_file_index", &source_file_index);
  pairs.SetBranchAddress("event_id_hi", &event_id_hi);
  pairs.SetBranchAddress("event_id_lo", &event_id_lo);
  pairs.SetBranchAddress("candidate_id_hi", &candidate_id_hi);
  pairs.SetBranchAddress("candidate_id_lo", &candidate_id_lo);
  pairs.SetBranchAddress("jet_radius", &jet_radius);
  pairs.SetBranchAddress("jet_pt", &jet_pt);
  pairs.SetBranchAddress("jet_eta", &jet_eta);
  pairs.SetBranchAddress("delta_phi", &delta_phi);
  pairs.SetBranchAddress("xjgamma", &xjgamma);

  // This package-only example is deliberately raw/unweighted. In particular,
  // event_weight is not the complete Au+Au inclusive producer*stitch*centrality
  // weight. Use a separately certified composite payload for normalized output.
  TH1D *hist = new TH1D("h_xjgamma_raw", ";x_{J#gamma}=p_{T}^{jet}/p_{T}^{#gamma};raw selected pairs", 18, 0, 1.8);
  for (Long64_t i = 0; i < pairs.GetEntries(); ++i) {
    pairs.GetEntry(i);
    const auto found = leaders.find(
      {pairs.GetTreeNumber(), source_file_index, event_id_hi, event_id_lo}
    );
    if (found == leaders.end()) continue;
    if (candidate_id_hi != found->second.candidate_id_hi ||
        candidate_id_lo != found->second.candidate_id_lo) continue;
    // Always select the jet radius explicitly. This prevents an additive
    // multi-radius release from silently mixing R=0.3 and R=0.4 jets.
    if (std::abs(jet_radius - 0.4) > 1e-12) continue;
    if (jet_pt <= 5 || std::abs(jet_eta) >= 0.7) continue;
    if (delta_phi <= 7.0 * TMath::Pi() / 8.0) continue;
    hist->Fill(xjgamma);
  }

  printf("selected event-leading tight, R=0.4-isolated photons: %zu\n", leaders.size());
  printf("RAW/UNWEIGHTED diagnostic: not a normalized Au+Au inclusive result\n");
  TCanvas *canvas = new TCanvas("c_xjgamma", "xJgamma quick start", 800, 650);
  hist->Draw("E");
  canvas->SaveAs("quick_start_xjgamma.png");
}
