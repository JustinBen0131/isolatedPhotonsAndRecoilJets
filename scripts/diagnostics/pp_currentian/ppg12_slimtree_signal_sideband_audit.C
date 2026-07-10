// Fast foreground diagnostic for THE-76 PPG12 Fig.29 leakage parity.
//
// Reads PPG12 bdt_split slimtrees and reproduces the signal-cluster
// tight/non-tight/iso/noniso classification used by
// efficiencytool/RecoEffCalculator_TTreeReader.C.  This is a diagnostic only:
// it does not write PPG12 outputs and is intended for SDCC scratch execution.

#include <TFile.h>
#include <TH2.h>
#include <TRandom3.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
constexpr int kMaxClusters = 20000;
constexpr int kMaxParticles = 20000;
constexpr double kRecoBins[] = {10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36};
constexpr int kNBins = static_cast<int>(sizeof(kRecoBins) / sizeof(kRecoBins[0])) - 1;

struct Counts
{
  double all = 0;
  double tight = 0;
  double A = 0;
  double B = 0;
  double C = 0;
  double D = 0;
  double common = 0;
  double nt = 0;
  double neither = 0;
  double signal_cluster = 0;
  double bdt_tight = 0;
  double bdt_nt_window = 0;
  double nt_shape = 0;
  double nt_fail_any = 0;
  double fail_weta = 0;
  double fail_prob = 0;
  double fail_bdt = 0;
  double fail_weta_only = 0;
  double fail_bdt_only = 0;
  double fail_prob_only = 0;
  double fail_weta_bdt = 0;
  double fail_other_combo = 0;
  double score_sum = 0;
  double score_sumw = 0;
};

int find_bin(const double x)
{
  for (int i = 0; i < kNBins; ++i)
  {
    if (x > kRecoBins[i] && x < kRecoBins[i + 1]) return i;
  }
  return -1;
}

bool open_interval(const double x, const double lo, const double hi)
{
  return std::isfinite(x) && x > lo && x < hi;
}

bool tower_masked(TH2* mask, const float ieta_value, const float iphi_value)
{
  if (!mask) return false;
  const int ieta = static_cast<int>(ieta_value);
  const int iphi = static_cast<int>(iphi_value);
  if (ieta < 0 || ieta >= mask->GetNbinsX()) return false;
  if (iphi < 0 || iphi >= mask->GetNbinsY()) return false;
  return mask->GetBinContent(ieta + 1, iphi + 1) > 0;
}

double safe_div(const double n, const double d)
{
  return d != 0 ? n / d : std::numeric_limits<double>::quiet_NaN();
}

void print_counts(const std::string& label, const std::vector<Counts>& c)
{
  std::cout << "\n" << label << "\n";
  std::cout << "bin,A,B,C,D,all,tight,common,nt,neither,B/A,C/A,D/A,tight/all,nt/all,"
            << "signal_cluster,common/signal,nt/common,neither/common,bdt_tight/common,bdt_nt_window/common,"
            << "nt_shape/common,nt_fail_any/common,fail_weta/common,fail_prob/common,fail_bdt/common,"
            << "fail_weta_only/common,fail_bdt_only/common,fail_prob_only/common,fail_weta_bdt/common,"
            << "fail_other_combo/common,mean_bdt\n";
  for (int i = 0; i < kNBins; ++i)
  {
    const auto& v = c[i];
    std::cout << std::fixed << std::setprecision(0)
              << kRecoBins[i] << "-" << kRecoBins[i + 1] << ","
              << std::setprecision(6)
              << v.A << "," << v.B << "," << v.C << "," << v.D << ","
              << v.all << "," << v.tight << "," << v.common << "," << v.nt << "," << v.neither << ","
              << safe_div(v.B, v.A) << "," << safe_div(v.C, v.A) << "," << safe_div(v.D, v.A) << ","
              << safe_div(v.tight, v.all) << "," << safe_div(v.nt, v.all) << ","
              << v.signal_cluster << ","
              << safe_div(v.common, v.signal_cluster) << ","
              << safe_div(v.nt, v.common) << ","
              << safe_div(v.neither, v.common) << ","
              << safe_div(v.bdt_tight, v.common) << ","
              << safe_div(v.bdt_nt_window, v.common) << ","
              << safe_div(v.nt_shape, v.common) << ","
              << safe_div(v.nt_fail_any, v.common) << ","
              << safe_div(v.fail_weta, v.common) << ","
              << safe_div(v.fail_prob, v.common) << ","
              << safe_div(v.fail_bdt, v.common) << ","
              << safe_div(v.fail_weta_only, v.common) << ","
              << safe_div(v.fail_bdt_only, v.common) << ","
              << safe_div(v.fail_prob_only, v.common) << ","
              << safe_div(v.fail_weta_bdt, v.common) << ","
              << safe_div(v.fail_other_combo, v.common) << ","
              << safe_div(v.score_sum, v.score_sumw)
              << "\n";
  }
}
} // namespace

void ppg12_slimtree_signal_sideband_audit(long long max_entries = -1);

void the76_ppg12_slimtree_signal_sideband_audit(long long max_entries = -1)
{
  ppg12_slimtree_signal_sideband_audit(max_entries);
}

void ppg12_slimtree_signal_sideband_audit(long long max_entries)
{
  const std::vector<std::string> files = {
      "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon5/bdt_split.root",
      "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon10/bdt_split.root",
      "/sphenix/user/shuhangli/ppg12/FunWithxgboost/photon20/bdt_split.root"};
  const std::vector<double> weights = {
      146359.3 / 130.4461,
      6944.675 / 130.4461,
      1.0};

  TFile mask_file("/sphenix/user/shuhangli/ppg12/efficiencytool/tower_masks_bdt_nom.root");
  TH2* mask = dynamic_cast<TH2*>(mask_file.Get("mask_phisymm_tight"));

  std::vector<Counts> total(kNBins);
  TRandom3 rand(0);

  for (std::size_t ifile = 0; ifile < files.size(); ++ifile)
  {
    TFile input(files[ifile].c_str());
    auto* tree = dynamic_cast<TTree*>(input.Get("slimtree"));
    if (!tree)
    {
      std::cerr << "missing slimtree in " << files[ifile] << "\n";
      continue;
    }

    float vertexz = 0;
    int nparticles = 0;
    float particle_Pt[kMaxParticles];
    float particle_Eta[kMaxParticles];
    int particle_pid[kMaxParticles];
    int particle_trkid[kMaxParticles];
    int particle_photonclass[kMaxParticles];
    int particle_converted[kMaxParticles];
    float particle_truth_iso_03[kMaxParticles];

    int ncluster = 0;
    float cluster_Et[kMaxClusters];
    float cluster_Eta[kMaxClusters];
    float cluster_prob[kMaxClusters];
    float cluster_iso_topo_04[kMaxClusters];
    float cluster_et1[kMaxClusters];
    float cluster_et2[kMaxClusters];
    float cluster_et3[kMaxClusters];
    float cluster_et4[kMaxClusters];
    float cluster_weta_cogx[kMaxClusters];
    float cluster_wphi_cogx[kMaxClusters];
    float cluster_ietacent[kMaxClusters];
    float cluster_iphicent[kMaxClusters];
    float cluster_e11[kMaxClusters];
    float cluster_e33[kMaxClusters];
    float cluster_e32[kMaxClusters];
    float cluster_e35[kMaxClusters];
    int cluster_truthtrkID[kMaxClusters];
    float cluster_npb_score[kMaxClusters];
    float bdt_base_E[kMaxClusters];
    float bdt_base_v3E[kMaxClusters];

    tree->SetBranchAddress("vertexz", &vertexz);
    tree->SetBranchAddress("nparticles", &nparticles);
    tree->SetBranchAddress("particle_Pt", particle_Pt);
    tree->SetBranchAddress("particle_Eta", particle_Eta);
    tree->SetBranchAddress("particle_pid", particle_pid);
    tree->SetBranchAddress("particle_trkid", particle_trkid);
    tree->SetBranchAddress("particle_photonclass", particle_photonclass);
    tree->SetBranchAddress("particle_converted", particle_converted);
    tree->SetBranchAddress("particle_truth_iso_03", particle_truth_iso_03);

    tree->SetBranchAddress("ncluster_CLUSTERINFO_CEMC", &ncluster);
    tree->SetBranchAddress("cluster_Et_CLUSTERINFO_CEMC", cluster_Et);
    tree->SetBranchAddress("cluster_Eta_CLUSTERINFO_CEMC", cluster_Eta);
    tree->SetBranchAddress("cluster_prob_CLUSTERINFO_CEMC", cluster_prob);
    tree->SetBranchAddress("cluster_iso_topo_04_CLUSTERINFO_CEMC", cluster_iso_topo_04);
    tree->SetBranchAddress("cluster_et1_CLUSTERINFO_CEMC", cluster_et1);
    tree->SetBranchAddress("cluster_et2_CLUSTERINFO_CEMC", cluster_et2);
    tree->SetBranchAddress("cluster_et3_CLUSTERINFO_CEMC", cluster_et3);
    tree->SetBranchAddress("cluster_et4_CLUSTERINFO_CEMC", cluster_et4);
    tree->SetBranchAddress("cluster_weta_cogx_CLUSTERINFO_CEMC", cluster_weta_cogx);
    tree->SetBranchAddress("cluster_wphi_cogx_CLUSTERINFO_CEMC", cluster_wphi_cogx);
    tree->SetBranchAddress("cluster_ietacent_CLUSTERINFO_CEMC", cluster_ietacent);
    tree->SetBranchAddress("cluster_iphicent_CLUSTERINFO_CEMC", cluster_iphicent);
    tree->SetBranchAddress("cluster_e11_CLUSTERINFO_CEMC", cluster_e11);
    tree->SetBranchAddress("cluster_e33_CLUSTERINFO_CEMC", cluster_e33);
    tree->SetBranchAddress("cluster_e32_CLUSTERINFO_CEMC", cluster_e32);
    tree->SetBranchAddress("cluster_e35_CLUSTERINFO_CEMC", cluster_e35);
    tree->SetBranchAddress("cluster_truthtrkID_CLUSTERINFO_CEMC", cluster_truthtrkID);
    tree->SetBranchAddress("cluster_npb_score_CLUSTERINFO_CEMC", cluster_npb_score);
    tree->SetBranchAddress("cluster_bdt_CLUSTERINFO_CEMC_base_E", bdt_base_E);
    tree->SetBranchAddress("cluster_bdt_CLUSTERINFO_CEMC_base_v3E", bdt_base_v3E);

    const auto entries = tree->GetEntries();
    const auto stop = (max_entries > 0) ? std::min<long long>(entries, max_entries) : entries;
    std::vector<Counts> sample(kNBins);

    for (long long ie = 0; ie < stop; ++ie)
    {
      tree->GetEntry(ie);
      if (std::abs(vertexz) > 60.0) continue;

      std::map<int, int> particle_by_track;
      std::map<int, bool> signal_particle;
      for (int ip = 0; ip < nparticles && ip < kMaxParticles; ++ip)
      {
        particle_by_track[particle_trkid[ip]] = ip;
        if (particle_pid[ip] == 22 &&
            particle_photonclass[ip] < 3 &&
            particle_truth_iso_03[ip] < 4.0)
        {
          signal_particle[ip] = true;
        }
      }

      std::vector<double> smeared_et(std::min(ncluster, kMaxClusters), 0.0);
      for (int ic = 0; ic < ncluster && ic < kMaxClusters; ++ic)
      {
        if (tower_masked(mask, cluster_ietacent[ic], cluster_iphicent[ic])) continue;
        smeared_et[ic] = cluster_Et[ic] * rand.Gaus(1.0, 0.04);
      }

      for (int ic = 0; ic < ncluster && ic < kMaxClusters; ++ic)
      {
        if (tower_masked(mask, cluster_ietacent[ic], cluster_iphicent[ic])) continue;
        const double et = smeared_et[ic];
        if (et < 5.0) continue;

        auto pit = particle_by_track.find(cluster_truthtrkID[ic]);
        if (pit == particle_by_track.end()) continue;
        const int ip = pit->second;
        if (!signal_particle.count(ip)) continue;

        const double e11e33 = cluster_e33[ic] > 0 ? cluster_e11[ic] / cluster_e33[ic] : 0.0;
        const double e32e35 = cluster_e35[ic] > 0 ? cluster_e32[ic] / cluster_e35[ic] : 0.0;
        const double wr = cluster_weta_cogx[ic] != 0.0 ? cluster_wphi_cogx[ic] / cluster_weta_cogx[ic] : -999.0;
        const double bdt = (et >= 8.0 && et < 35.0) ? bdt_base_v3E[ic] : bdt_base_E[ic];
        const bool common =
            open_interval(cluster_prob[ic], 0.0, 1.0) &&
            open_interval(e11e33, 0.0, 0.98) &&
            wr > 0.0 &&
            cluster_weta_cogx[ic] < 2.0 &&
            cluster_npb_score[ic] > 0.5;

        const double tight_bdt_min = 0.815625 - 0.0015625 * et;
        const double nt_bdt_min = 0.7333333333333333 - 0.01333333333333333 * et;
        const double nt_bdt_max = 0.684375 + 0.0015625 * et;

        const bool tightProb = open_interval(cluster_prob[ic], 0.0, 1.0);
        const bool tightWeta = open_interval(cluster_weta_cogx[ic], 0.0, 1.0);
        const bool tightWphi = open_interval(cluster_wphi_cogx[ic], 0.0, 1.0);
        const bool tightEt1 = open_interval(cluster_et1[ic], 0.5, 1.0);
        const bool tightEt2 = open_interval(cluster_et2[ic], 0.0, 1.0);
        const bool tightEt3 = open_interval(cluster_et3[ic], 0.0, 1.0);
        const bool tightEt4 = open_interval(cluster_et4[ic], 0.0, 1.0);
        const bool tightE11E33 = open_interval(e11e33, 0.0, 1.0);
        const bool tightE32E35 = open_interval(e32e35, 0.8, 1.0);
        const bool tightBDT = bdt > tight_bdt_min && bdt < 1.0;
        const bool tight = common && tightProb && tightWeta && tightWphi && tightEt1 && tightEt2 &&
                           tightEt3 && tightEt4 && tightE11E33 && tightE32E35 && tightBDT;

        const bool ntShape =
            common &&
            open_interval(cluster_prob[ic], 0.0, 1.0) &&
            open_interval(cluster_weta_cogx[ic], 0.0, 1.0) &&
            open_interval(cluster_wphi_cogx[ic], 0.0, 1.0) &&
            open_interval(cluster_et1[ic], 0.6, 1.0) &&
            open_interval(cluster_et4[ic], 0.0, 1.0) &&
            open_interval(e11e33, 0.0, 1.0) &&
            open_interval(e32e35, 0.8, 1.0);
        int nfail = 0;
        if (!tightWeta) ++nfail;
        if (!tightProb) ++nfail;
        if (!tightBDT) ++nfail;
        const bool nontight = ntShape && bdt > nt_bdt_min && bdt < nt_bdt_max && nfail > 0;

        const double eiso = cluster_iso_topo_04[ic] * 1.2 + 0.1;
        const double iso_max = 0.490 + 0.037 * et;
        const bool iso = eiso > -20.0 && eiso < iso_max;
        const bool noniso = eiso > iso_max + 0.8 && eiso < 20.0;
        const int ib = find_bin(et);
        const double w = weights[ifile];
        if (ib >= 0)
        {
          sample[ib].signal_cluster += w;
          total[ib].signal_cluster += w;
          sample[ib].common += common ? w : 0.0;
          total[ib].common += common ? w : 0.0;
          if (common)
          {
            sample[ib].score_sum += bdt * w;
            sample[ib].score_sumw += w;
            total[ib].score_sum += bdt * w;
            total[ib].score_sumw += w;
            if (tightBDT)
            {
              sample[ib].bdt_tight += w;
              total[ib].bdt_tight += w;
            }
            if (bdt > nt_bdt_min && bdt < nt_bdt_max)
            {
              sample[ib].bdt_nt_window += w;
              total[ib].bdt_nt_window += w;
            }
            if (ntShape)
            {
              sample[ib].nt_shape += w;
              total[ib].nt_shape += w;
            }
            if (nfail > 0)
            {
              sample[ib].nt_fail_any += w;
              total[ib].nt_fail_any += w;
            }
            const bool failWeta = !tightWeta;
            const bool failProb = !tightProb;
            const bool failBdt = !tightBDT;
            if (failWeta)
            {
              sample[ib].fail_weta += w;
              total[ib].fail_weta += w;
            }
            if (failProb)
            {
              sample[ib].fail_prob += w;
              total[ib].fail_prob += w;
            }
            if (failBdt)
            {
              sample[ib].fail_bdt += w;
              total[ib].fail_bdt += w;
            }
            const int failMask =
                (failWeta ? 1 : 0) |
                (failProb ? 2 : 0) |
                (failBdt ? 4 : 0);
            if (failMask == 1)
            {
              sample[ib].fail_weta_only += w;
              total[ib].fail_weta_only += w;
            }
            else if (failMask == 2)
            {
              sample[ib].fail_prob_only += w;
              total[ib].fail_prob_only += w;
            }
            else if (failMask == 4)
            {
              sample[ib].fail_bdt_only += w;
              total[ib].fail_bdt_only += w;
            }
            else if (failMask == 5)
            {
              sample[ib].fail_weta_bdt += w;
              total[ib].fail_weta_bdt += w;
            }
            else if (failMask != 0)
            {
              sample[ib].fail_other_combo += w;
              total[ib].fail_other_combo += w;
            }
          }
          if (particle_Pt[ip] > 8.0 && particle_Pt[ip] < 45.0 && et > 10.0 && et < 36.0)
          {
            sample[ib].all += w;
            total[ib].all += w;
          }
          if (tight)
          {
            if (particle_Pt[ip] > 8.0 && particle_Pt[ip] < 45.0 && et > 10.0 && et < 36.0)
            {
              sample[ib].tight += w;
              total[ib].tight += w;
              if (iso)
              {
                sample[ib].A += w;
                total[ib].A += w;
              }
            }
            if (noniso)
            {
              sample[ib].B += w;
              total[ib].B += w;
            }
          }
          else if (nontight)
          {
            sample[ib].nt += w;
            total[ib].nt += w;
            if (iso)
            {
              sample[ib].C += w;
              total[ib].C += w;
            }
            if (noniso)
            {
              sample[ib].D += w;
              total[ib].D += w;
            }
          }
          else if (common)
          {
            sample[ib].neither += w;
            total[ib].neither += w;
          }
        }
      }
    }
    std::cout << "processed " << stop << " entries from " << files[ifile] << "\n";
    print_counts(files[ifile], sample);
  }

  print_counts("weighted total", total);
}
