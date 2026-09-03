#include <TError.h>
#include <TChain.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <tuple>

struct RegionLeader {
  double photon_et = -1;
  int encounter_ordinal = -1;
};

static bool purity_add_files(TChain &chain, const char *path)
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

// Choose non_tight_definition explicitly:
//   "bounded" uses bdt_is_nontight;
//   "complement" uses bdt_is_not_tight.
// isolation_radius may be 0.3 or 0.4. Examples:
//   purity_regions("../data/auau/files.txt", "auau", "bounded");
//   purity_regions("../data/auau/files.txt", "auau", "complement");
void purity_regions(const char *product,
                    const char *system,
                    const char *non_tight_definition,
                    double isolation_radius = 0.4)
{
  const std::string collision_system = system ? system : "";
  if (collision_system != "pp" && collision_system != "auau") {
    Error("purity_regions", "system must be pp or auau");
    return;
  }
  const std::string definition = non_tight_definition ? non_tight_definition : "";
  if (definition != "bounded" && definition != "complement") {
    Error("purity_regions", "non_tight_definition must be bounded or complement");
    return;
  }
  const bool use_bounded_non_tight = definition == "bounded";
  if (isolation_radius != 0.3 && isolation_radius != 0.4) {
    Error("purity_regions", "isolation_radius must be 0.3 or 0.4");
    return;
  }

  TChain tree("photons");
  if (!purity_add_files(tree, product) || tree.GetEntries() < 0) {
    Error("purity_regions", "Cannot load photons from %s", product); return;
  }

  Int_t source_file_index = -1;
  ULong64_t event_id_hi = 0, event_id_lo = 0;
  Int_t encounter_ordinal = -1;
  double photon_et = 0, photon_eta = 0, vertex_z = 0, centrality = 0;
  Int_t tight = 0, bounded_nontight = 0, complement_nontight = 0;
  double isolation = 0, isolated_threshold = 0, nonisolated_threshold = 0;

  const std::string radius = isolation_radius == 0.3 ? "r03" : "r04";

  tree.SetBranchAddress("source_file_index", &source_file_index);
  tree.SetBranchAddress("event_id_hi", &event_id_hi);
  tree.SetBranchAddress("event_id_lo", &event_id_lo);
  tree.SetBranchAddress("photon_encounter_ordinal", &encounter_ordinal);
  tree.SetBranchAddress("photon_et", &photon_et);
  tree.SetBranchAddress("photon_eta", &photon_eta);
  tree.SetBranchAddress("vertex_z", &vertex_z);
  tree.SetBranchAddress("centrality", &centrality);
  tree.SetBranchAddress("bdt_is_tight", &tight);
  tree.SetBranchAddress("bdt_is_nontight", &bounded_nontight);
  tree.SetBranchAddress("bdt_is_not_tight", &complement_nontight);
  tree.SetBranchAddress(("iso_" + radius).c_str(), &isolation);
  tree.SetBranchAddress(("iso_" + radius + "_threshold").c_str(), &isolated_threshold);
  tree.SetBranchAddress(("iso_" + radius + "_nonisolated_threshold").c_str(), &nonisolated_threshold);

  using EventKey = std::tuple<Int_t, Int_t, ULong64_t, ULong64_t>;
  std::map<EventKey, std::array<RegionLeader, 4>> leaders;
  for (Long64_t i = 0; i < tree.GetEntries(); ++i) {
    tree.GetEntry(i);
    if (photon_et < 15 || photon_et >= 35) continue;
    if (std::abs(photon_eta) >= 0.7) continue;
    if (collision_system == "pp" && std::abs(vertex_z) >= 30) continue;
    if (collision_system == "auau" &&
        (std::abs(vertex_z) >= 10 || centrality < 0 || centrality >= 20)) continue;

    const bool non_tight = use_bounded_non_tight ? bounded_nontight : complement_nontight;
    const bool isolated = isolation < isolated_threshold;
    const bool nonisolated = isolation > nonisolated_threshold;

    int region = -1;
    if (tight && isolated) region = 0;
    if (tight && nonisolated) region = 1;
    if (non_tight && isolated) region = 2;
    if (non_tight && nonisolated) region = 3;
    if (region < 0) continue;

    RegionLeader &leader = leaders[
      {tree.GetTreeNumber(), source_file_index, event_id_hi, event_id_lo}
    ][region];
    const bool precedes = leader.photon_et < 0 || photon_et > leader.photon_et ||
      (photon_et == leader.photon_et && encounter_ordinal < leader.encounter_ordinal);
    if (precedes) leader = {photon_et, encounter_ordinal};
  }

  Long64_t counts[4] = {0, 0, 0, 0};
  for (const auto &event : leaders)
    for (int region = 0; region < 4; ++region)
      if (event.second[region].photon_et >= 0) counts[region] += 1;

  printf("system: %s; 15 <= photon pT < 35 GeV; |photon eta| < 0.7; ",
         collision_system.c_str());
  if (collision_system == "pp") printf("|vertex z| < 30 cm\n");
  else printf("|vertex z| < 10 cm; 0 <= centrality < 20%%\n");
  printf("non-tight definition: %s; isolation radius: %.1f\n",
         use_bounded_non_tight ? "bounded score band" : "full tight-score complement",
         isolation_radius);
  printf("raw/unweighted event-leading A=%lld B=%lld C=%lld D=%lld\n",
         counts[0], counts[1], counts[2], counts[3]);
  printf("Package-only diagnostic: not a normalized Au+Au inclusive result.\n");
}
