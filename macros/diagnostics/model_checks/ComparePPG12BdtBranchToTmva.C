#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TMVA/RBDT.hxx>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
float safe_ratio(float num, float den)
{
  return den > 0.0F ? num / den : 0.0F;
}

bool finite_all(const std::vector<float>& values)
{
  for (const float value : values)
  {
    if (!std::isfinite(value)) return false;
  }
  return true;
}
}  // namespace

void ComparePPG12BdtBranchToTmva(const char* sample = "photon20", long long maxEntries = 20000)
{
  const std::string treePath =
      std::string("/sphenix/user/shuhangli/ppg12/FunWithxgboost/") + sample + "/bdt_split.root";
  const std::string modelPath =
      "/sphenix/user/shuhangli/ppg12/FunWithxgboost/binned_models/model_base_v3E_split_single_tmva.root";
  const std::string node = "CLUSTERINFO_CEMC";

  TFile input(treePath.c_str(), "READ");
  if (input.IsZombie())
  {
    std::cerr << "ERROR opening " << treePath << std::endl;
    return;
  }
  auto* tree = dynamic_cast<TTree*>(input.Get("slimtree"));
  if (!tree)
  {
    std::cerr << "ERROR: missing slimtree in " << treePath << std::endl;
    return;
  }

  TMVA::Experimental::RBDT model("myBDT", modelPath);

  TTreeReader reader(tree);
  TTreeReaderValue<float> vertexz(reader, "vertexz");
  TTreeReaderArray<float> cluster_Et(reader, ("cluster_Et_" + node).c_str());
  TTreeReaderArray<float> cluster_Eta(reader, ("cluster_Eta_" + node).c_str());
  TTreeReaderArray<float> cluster_weta_cogx(reader, ("cluster_weta_cogx_" + node).c_str());
  TTreeReaderArray<float> cluster_wphi_cogx(reader, ("cluster_wphi_cogx_" + node).c_str());
  TTreeReaderArray<float> cluster_e11(reader, ("cluster_e11_" + node).c_str());
  TTreeReaderArray<float> cluster_e33(reader, ("cluster_e33_" + node).c_str());
  TTreeReaderArray<float> cluster_e32(reader, ("cluster_e32_" + node).c_str());
  TTreeReaderArray<float> cluster_e35(reader, ("cluster_e35_" + node).c_str());
  TTreeReaderArray<float> cluster_et1(reader, ("cluster_et1_" + node).c_str());
  TTreeReaderArray<float> cluster_et2(reader, ("cluster_et2_" + node).c_str());
  TTreeReaderArray<float> cluster_et3(reader, ("cluster_et3_" + node).c_str());
  TTreeReaderArray<float> cluster_et4(reader, ("cluster_et4_" + node).c_str());
  TTreeReaderArray<float> cluster_bdt(reader, ("cluster_bdt_" + node + "_base_v3E").c_str());

  long long entries = 0;
  long long clusters = 0;
  long long selected = 0;
  long long badFinite = 0;
  long long over1e6 = 0;
  long long over1e5 = 0;
  long long over1e4 = 0;
  double sumDiff = 0.0;
  double sumAbsDiff = 0.0;
  double sumSqDiff = 0.0;
  double maxAbsDiff = -1.0;
  long long maxEntry = -1;
  int maxCluster = -1;
  float maxStored = std::numeric_limits<float>::quiet_NaN();
  float maxComputed = std::numeric_limits<float>::quiet_NaN();

  std::cout << std::fixed << std::setprecision(9);
  while (reader.Next())
  {
    if (maxEntries >= 0 && entries >= maxEntries) break;
    const long long entry = entries;
    ++entries;
    const auto nCluster = cluster_Et.GetSize();
    clusters += nCluster;
    for (int ic = 0; ic < nCluster; ++ic)
    {
      const float et = cluster_Et[ic];
      if (!(et >= 8.0F && et < 35.0F)) continue;

      const std::vector<float> x = {
          et,
          cluster_weta_cogx[ic],
          cluster_wphi_cogx[ic],
          *vertexz,
          cluster_Eta[ic],
          safe_ratio(cluster_e11[ic], cluster_e33[ic]),
          cluster_et1[ic],
          cluster_et2[ic],
          cluster_et3[ic],
          cluster_et4[ic],
          safe_ratio(cluster_e32[ic], cluster_e35[ic])
      };
      if (!finite_all(x) || !std::isfinite(cluster_bdt[ic]))
      {
        ++badFinite;
        continue;
      }

      const float computed = model.Compute(x)[0];
      const float stored = cluster_bdt[ic];
      const double diff = static_cast<double>(computed) - stored;
      const double absDiff = std::fabs(diff);
      ++selected;
      sumDiff += diff;
      sumAbsDiff += absDiff;
      sumSqDiff += diff * diff;
      if (absDiff > 1e-6) ++over1e6;
      if (absDiff > 1e-5) ++over1e5;
      if (absDiff > 1e-4) ++over1e4;
      if (absDiff > maxAbsDiff)
      {
        maxAbsDiff = absDiff;
        maxEntry = entry;
        maxCluster = ic;
        maxStored = stored;
        maxComputed = computed;
      }
      if (absDiff > 1e-5 && over1e5 <= 8)
      {
        std::cout << "BAD_ROW"
                  << " sample=" << sample
                  << " entry=" << entry
                  << " cluster=" << ic
                  << " et=" << et
                  << " eta=" << cluster_Eta[ic]
                  << " stored=" << stored
                  << " computed=" << computed
                  << " diff=" << diff
                  << " weta=" << cluster_weta_cogx[ic]
                  << " wphi=" << cluster_wphi_cogx[ic]
                  << " vtx=" << *vertexz
                  << " e11e33=" << x[5]
                  << " e32e35=" << x[10]
                  << std::endl;
      }
    }
  }

  const double meanDiff = selected > 0 ? sumDiff / selected : 0.0;
  const double meanAbsDiff = selected > 0 ? sumAbsDiff / selected : 0.0;
  const double rmsDiff = selected > 0 ? std::sqrt(sumSqDiff / selected) : 0.0;
  std::cout << "SUMMARY"
            << " sample=" << sample
            << " entries=" << entries
            << " clusters=" << clusters
            << " selected_8to35=" << selected
            << " bad_finite=" << badFinite
            << " mean_diff=" << meanDiff
            << " mean_abs_diff=" << meanAbsDiff
            << " rms_diff=" << rmsDiff
            << " max_abs_diff=" << maxAbsDiff
            << " over_1e-6=" << over1e6
            << " over_1e-5=" << over1e5
            << " over_1e-4=" << over1e4
            << " max_entry=" << maxEntry
            << " max_cluster=" << maxCluster
            << " max_stored=" << maxStored
            << " max_computed=" << maxComputed
            << std::endl;
}
