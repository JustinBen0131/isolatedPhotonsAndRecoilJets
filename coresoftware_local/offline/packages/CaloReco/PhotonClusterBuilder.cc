#include "PhotonClusterBuilder.h"

#include <calobase/PhotonClusterv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <calobase/TowerInfoDefs.h>

// Tower stuff
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfo.h>

// for the vertex
#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <TMVA/RBDT.hxx>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>   // NEW: for std::setw / std::setprecision in debug tables
#include <set>
#include <stdexcept>

namespace
{
  // Helper function to shift tower indices for wrapping in phi
  void shift_tower_index(int& ieta, int& iphi, int etadiv, int phidiv)
  {
    while (iphi < 0)
    {
      iphi += phidiv;
    }
    while (iphi >= phidiv)
    {
      iphi -= phidiv;
    }
    if (ieta < 0 || ieta >= etadiv)
    {
      ieta = -1;  // invalid
    }
  }
}  // namespace

PhotonClusterBuilder::PhotonClusterBuilder(const std::string& name)
  : SubsysReco(name)
{
}

int PhotonClusterBuilder::InitRun(PHCompositeNode* topNode)
{
  // BDT
  if (m_do_bdt)
  {
    m_bdt = std::make_unique<TMVA::Experimental::RBDT>("myBDT", m_bdt_model_file);
  }

  // locate input raw cluster container
  m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
  if (!m_rawclusters)
  {
    std::cerr << Name() << ": could not find RawClusterContainer node '" << m_input_cluster_node << "'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_emc_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  if (!m_emc_tower_container)
  {
    std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_CEMC'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  if (!m_geomEM)
    {
      m_geomEM = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC_DETAILED");
    }

    if (!m_geomEM)
    {
      std::cerr << Name()
                << ": could not find RawTowerGeomContainer node 'TOWERGEOM_CEMC' "
                << "(or fallback 'TOWERGEOM_CEMC_DETAILED')"
                << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ihcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  if (!m_ihcal_tower_container)
  {
    std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_HCALIN'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  if (!m_geomIH)
  {
    std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALIN'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ohcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if (!m_ohcal_tower_container)
  {
    std::cerr << Name() << ": could not find TowerInfoContainer node 'TOWERINFO_CALIB_HCALOUT'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  if (!m_geomOH)
  {
    std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALOUT'" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonClusterBuilder::CreateNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    throw std::runtime_error("PhotonClusterBuilder: DST node not found");
  }
  m_photon_container = findNode::getClass<RawClusterContainer>(dstNode, m_output_photon_node);
  if (!m_photon_container)
  {
    m_photon_container = new RawClusterContainer();
    auto* photonNode = new PHIODataNode<PHObject>(m_photon_container, m_output_photon_node, "PHObject");
    dstNode->addNode(photonNode);
  }
}

int PhotonClusterBuilder::process_event(PHCompositeNode* topNode)
{
    if (!m_rawclusters)
    {
      m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
      if (!m_rawclusters)
      {
        std::cerr << Name() << ": missing RawClusterContainer '" << m_input_cluster_node << "'" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }

    // Always clear output container so skipped events cannot leave stale photons
    if (m_photon_container)
    {
      m_photon_container->Reset();
    }

    // ------------------------------------------------------------------
    // Vertex selection (RECO vertex for RECO objects)
    //   Use MBD vertex (DATA + SIM). Optional fallback: GlobalVertexMap (still reco).
    //   Never overwrite reco objects with TRUTH vertex.
    // ------------------------------------------------------------------
    static unsigned long long s_evt = 0;
    ++s_evt;

    const float vzCutForInfo = m_vz_cut_cm;

    m_vertex = std::numeric_limits<float>::quiet_NaN();
    const char* vtx_source = "NONE";

    // Snapshot what the maps contain (so we can print on skip)
    float mbd_z = std::numeric_limits<float>::quiet_NaN();
    float gv_z  = std::numeric_limits<float>::quiet_NaN();
    size_t mbd_n = 0;
    size_t gv_n  = 0;

    if (auto* mbdmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap"))
    {
      mbd_n = mbdmap->size();
      if (!mbdmap->empty() && mbdmap->begin()->second)
      {
        mbd_z = mbdmap->begin()->second->get_z();
      }
    }

    if (auto* gvmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap"))
    {
      gv_n = gvmap->size();
      if (!gvmap->empty() && gvmap->begin()->second)
      {
        gv_z = gvmap->begin()->second->get_z();
      }
    }

    // 1) MBD vertex (preferred, but only if finite)
    if (std::isfinite(mbd_z))
    {
      m_vertex = mbd_z;
      vtx_source = "MBD";
    }

    // 2) Optional fallback: GlobalVertexMap (still reco, only if finite)
    if (!std::isfinite(m_vertex) && std::isfinite(gv_z))
    {
      m_vertex = gv_z;
      vtx_source = "GlobalVertex";
    }

    // If still invalid, skip
    if (!std::isfinite(m_vertex))
    {
      if (Verbosity() >= 1)
      {
        std::cout << Name()
                  << ": [evt=" << s_evt << "] SKIP (no finite reco vertex)"
                  << " | MbdVertexMap n=" << mbd_n << " z=" << (std::isfinite(mbd_z) ? std::to_string(mbd_z) : std::string("NaN/NA"))
                  << " | GlobalVertexMap n=" << gv_n << " z=" << (std::isfinite(gv_z) ? std::to_string(gv_z) : std::string("NaN/NA"))
                  << std::endl;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }

    const bool passVz = (!m_use_vz_cut) || (std::fabs(m_vertex) < vzCutForInfo);

    if (Verbosity() >= 2)
    {
      std::cout << Name()
                << ": [evt=" << s_evt << "] using vertex_z=" << m_vertex
                << " (source=" << vtx_source << ")";
      if (m_use_vz_cut)
      {
        std::cout << " | inside |vz|<" << vzCutForInfo << "? " << (passVz ? "YES" : "NO");
      }
      else
      {
        std::cout << " | vz cut DISABLED";
      }
      std::cout << std::endl;
    }

    // If vertex is outside your analysis vz window, skip building photons for this event.
    if (m_use_vz_cut && !passVz)
    {
      if (Verbosity() >= 1)
      {
        std::cout << Name()
                  << ": [evt=" << s_evt << "] SKIP (|vz| cut)"
                  << " | vz=" << m_vertex << " | cut=" << vzCutForInfo
                  << " | source=" << vtx_source
                  << std::endl;
      }
      return Fun4AllReturnCodes::EVENT_OK;
    }

  // iterate over clusters via map to have access to keys if needed
  const auto& rcmap = m_rawclusters->getClustersMap();

  size_t nClusters = rcmap.size();
  size_t nPassET = 0;
  size_t nBuilt = 0;
  size_t nMeanTimeDefault = 0;

  for (const auto& kv : rcmap)
  {
    RawCluster* rc = kv.second;
    if (!rc)
    {
      continue;
    }

    CLHEP::Hep3Vector vertex_vec(0, 0, m_vertex);



    float eta = RawClusterUtility::GetPseudorapidity(*rc, vertex_vec);
    float phi = RawClusterUtility::GetAzimuthAngle(*rc, vertex_vec);
    float E = rc->get_energy();
    float ET = E / std::cosh(eta);
      if (ET < m_min_cluster_et)
      {
        continue;
      }
      ++nPassET;

        PhotonClusterv1* photon = new PhotonClusterv1(*rc);

      // ------------------------------------------------------------------
      // Persist kinematics used by PhotonClusterBuilder so downstream
      // (e.g. RecoilJets) can read consistent eta/phi/pt and vertex.
      //
      // NOTE:
      //   - We compute eta/phi using the chosen vertex_z (truth for SIM, MBD for DATA)
      //   - ET = E / cosh(eta)
      //   - For photons, pT == ET (massless)
      // ------------------------------------------------------------------
      photon->set_shower_shape_parameter("vertex_z",    m_vertex);
      photon->set_shower_shape_parameter("cluster_eta", eta);
      photon->set_shower_shape_parameter("cluster_phi", phi);
      photon->set_shower_shape_parameter("cluster_et",  ET);
      photon->set_shower_shape_parameter("cluster_pt",  ET);

      calculate_shower_shapes(rc, photon, eta, phi);

      const float mt = photon->get_shower_shape_parameter("mean_time");
      if (mt <= -998.5f) ++nMeanTimeDefault;

      //this is defensive coding, if do bdt is set false the bdt object should be nullptr
      //and this method will simply pass
      if (m_do_bdt)
      {
        calculate_bdt_score(photon);
      }

      ++nBuilt;
      m_photon_container->AddCluster(photon);
  }

  if (Verbosity() >= 2)
  {
    std::cout << Name()
              << ": [evt=" << s_evt << "] SUMMARY"
              << " | clusters=" << nClusters
              << " | passET(ET>" << m_min_cluster_et << ")=" << nPassET
              << " | photonsBuilt=" << nBuilt
              << " | mean_time=-999 count=" << nMeanTimeDefault
              << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void PhotonClusterBuilder::calculate_bdt_score(PhotonClusterv1* photon)
{
  if (!m_bdt)
  {
    return;
  }

  std::vector<float> x;
  for (const auto& feature : m_bdt_feature_list)
  {
    if (feature == "vertex_z")
    {
      x.push_back(m_vertex);
    }
    else if (feature == "ET")
    {
      float E = photon->get_energy();
      float ET = E / std::cosh(photon->get_shower_shape_parameter("cluster_eta"));
      x.push_back(ET);
    }
    else
    {
      const std::string delim = "_over_";
      const auto pos = feature.find(delim);
      if (pos != std::string::npos)
      {
        const std::string num = feature.substr(0, pos);
        const std::string den = feature.substr(pos + delim.size());
        const float numerator = photon->get_shower_shape_parameter(num);
        const float denominator = photon->get_shower_shape_parameter(den);
        x.push_back((denominator > 0) ? (numerator / denominator) : 0.0F);
      }
      else
      {
        x.push_back(photon->get_shower_shape_parameter(feature));
      }
    }
    //check if the thing we pushed back is NaN
    if (std::isnan(x.back()))
    {
      std::cerr << "PhotonClusterBuilder - feature name: " << feature << " is NaN" << std::endl;
    }
  }

  float bdt_score = -1;  // default value

  bdt_score = m_bdt->Compute(x)[0];

  photon->set_shower_shape_parameter("bdt_score", bdt_score);
}

void PhotonClusterBuilder::calculate_shower_shapes(RawCluster* rc, PhotonClusterv1* photon, float cluster_eta, float cluster_phi)
{
  std::vector<float> showershape = rc->get_shower_shapes(m_shape_min_tower_E);
  if (showershape.empty())
  {
    return;
  }

  std::pair<int, int> leadtowerindex = rc->get_lead_tower();
  int lead_ieta = leadtowerindex.first;
  int lead_iphi = leadtowerindex.second;

  float avg_eta = showershape[4] + 0.5F;
  float avg_phi = showershape[5] + 0.5F;


  int maxieta = std::floor(avg_eta);
  int maxiphi = std::floor(avg_phi);

  //if (maxieta < 3 || maxieta > 92)
  //{
  //  return;
  //}

  // for detamax, dphimax, nsaturated
  int detamax = 0;
  int dphimax = 0;
  int nsaturated = 0;
  float clusteravgtime = 0;
  float cluster_total_e = 0;
  const RawCluster::TowerMap& tower_map = rc->get_towermap();
  std::set<unsigned int> towers_in_cluster;
  for (auto tower_iter : tower_map)
  {
    RawTowerDefs::keytype tower_key = tower_iter.first;
    int ieta = RawTowerDefs::decode_index1(tower_key);
    int iphi = RawTowerDefs::decode_index2(tower_key);

    
    unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
    towers_in_cluster.insert(towerinfokey);
    TowerInfo* towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);
    if (towerinfo)
    {
      clusteravgtime += towerinfo->get_time() * towerinfo->get_energy();
      cluster_total_e += towerinfo->get_energy();
      if (towerinfo->get_isSaturated())
      {
        nsaturated++;
      }
    }
    
    int totalphibins = 256;
    auto dphiwrap = [totalphibins](int towerphi, int maxiphi_arg)
    {
      int idphi = towerphi - maxiphi_arg;
      if (idphi > totalphibins / 2)
      {
        idphi -= totalphibins;
      }
      if (idphi < -totalphibins / 2)
      {
        idphi += totalphibins;
      }
      return idphi;
    };

    int deta = ieta - lead_ieta;
    int dphi_val = dphiwrap(iphi, lead_iphi);

    detamax = std::max(std::abs(deta), detamax);
    dphimax = std::max(std::abs(dphi_val), dphimax);
  }

    if (cluster_total_e > 0)
    {
      clusteravgtime /= cluster_total_e;
    }
    else
    {
      if (Verbosity() >= 2)
      {
        std::cout << Name()
                  << ": cluster_total_e is 0, setting mean_time to -999"
                  << " | rc_energy=" << rc->get_energy()
                  << " | towers_in_towermap=" << tower_map.size()
                  << std::endl;

          // Diagnose WHY total_e stayed 0 (print a small sample)
          int nFound = 0;
          int nNonzero = 0;
          int nPrinted = 0;

          const float Epho = photon->get_energy();
          const float ETpho = Epho / std::cosh(cluster_eta);

          std::cout << "  context: E=" << Epho
                    << " ET=" << ETpho
                    << " eta=" << cluster_eta
                    << " phi=" << cluster_phi
                    << " lead(ieta,iphi)=(" << lead_ieta << "," << lead_iphi << ")"
                    << std::endl;

          for (auto tower_iter : tower_map)
          {
            RawTowerDefs::keytype tower_key = tower_iter.first;
            const float weight = tower_iter.second;

            int ieta = RawTowerDefs::decode_index1(tower_key);
            int iphi = RawTowerDefs::decode_index2(tower_key);

            unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
            TowerInfo* towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);

            if (towerinfo) nFound++;
            if (towerinfo && towerinfo->get_energy() > 0) nNonzero++;

            if (nPrinted < 6)
            {
              if (!towerinfo)
              {
                std::cout << "  tow[" << nPrinted << "] (ieta,iphi)=(" << ieta << "," << iphi << ")"
                          << " TI=MISSING"
                          << " raw_key=" << tower_key
                          << " ti_key=" << towerinfokey
                          << " weight=" << weight
                          << std::endl;
              }
              else
              {
                std::cout << "  tow[" << nPrinted << "] (ieta,iphi)=(" << ieta << "," << iphi << ")"
                          << " TI=OK"
                          << " isGood=" << (towerinfo->get_isGood() ? 1 : 0)
                          << " E=" << towerinfo->get_energy()
                          << " time=" << towerinfo->get_time()
                          << " raw_key=" << tower_key
                          << " ti_key=" << towerinfokey
                          << " weight=" << weight
                          << std::endl;
              }
              nPrinted++;
            }
          }

          std::cout << "  towerinfo summary: found=" << nFound
                    << " nonzeroE=" << nNonzero
                    << " (container size=" << m_emc_tower_container->size() << ")"
                    << std::endl;
      }

      clusteravgtime = -999.0f;
    }

  float E77[7][7] = {{0.0F}};
  int E77_ownership[7][7] = {{0}};

  for (int ieta = maxieta - 3; ieta < maxieta + 4; ieta++)
  {
    for (int iphi = maxiphi - 3; iphi < maxiphi + 4; iphi++)
    {
      //this is defensive coding, if ieta is out of range, set the energy to 0
      //even without this, the requirement for towerinfo object will take care of it
      if (ieta < 0 || ieta > 95)
      {
        E77[ieta - maxieta + 3][iphi - maxiphi + 3] = 0.0F;
        E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 0;
        continue;
      }

      int temp_ieta = ieta;
      int temp_iphi = iphi;
      shift_tower_index(temp_ieta, temp_iphi, 96, 256);
      if (temp_ieta < 0)
      {
        continue;
      }

      unsigned int towerinfokey = TowerInfoDefs::encode_emcal(temp_ieta, temp_iphi);

      if (towers_in_cluster.contains(towerinfokey))
      {
        E77_ownership[ieta - maxieta + 3][iphi - maxiphi + 3] = 1;
      }

      TowerInfo* towerinfo = m_emc_tower_container->get_tower_at_key(towerinfokey);
      if (towerinfo && towerinfo->get_isGood())
      {
        float energy = towerinfo->get_energy();
        if (energy > m_shape_min_tower_E)
        {
          E77[ieta - maxieta + 3][iphi - maxiphi + 3] = energy;
        }
      }
    }
  }

  float e11 = E77[3][3];
  float e33 = 0;
  float e55 = 0;
  float e77 = 0;
  float e13 = 0;
  float e15 = 0;
  float e17 = 0;
  float e31 = 0;
  float e51 = 0;
  float e71 = 0;
  float e35 = 0;
  float e37 = 0;
  float e53 = 0;
  float e73 = 0;
  float e57 = 0;
  float e75 = 0;
  float weta = 0;
  float wphi = 0;
  float weta_cog = 0;
  float wphi_cog = 0;
  float weta_cogx = 0;
  float wphi_cogx = 0;
  float Eetaphi = 0;
  float shift_eta = avg_eta - std::floor(avg_eta) - 0.5;
  float shift_phi = avg_phi - std::floor(avg_phi) - 0.5;
  float cog_eta = 3 + shift_eta;
  float cog_phi = 3 + shift_phi;

  float w32 = 0;
  float e32 = 0;
  float w52 = 0;
  float e52 = 0;
  float w72 = 0;
  float e72 = 0;
  float detacog = std::abs(maxieta - avg_eta);
  float dphicog = std::abs(maxiphi - avg_phi);
  float drad = std::sqrt(dphicog*dphicog + detacog*detacog);


  int signphi = (avg_phi - std::floor(avg_phi)) > 0.5 ? 1 : -1;

  for (int i = 0; i < 7; i++)
  {
    for (int j = 0; j < 7; j++)
    {
      int di = std::abs(i - 3);
      int dj = std::abs(j - 3);
      float di_float = i - cog_eta;
      float dj_float = j - cog_phi;

      if (E77_ownership[i][j] == 1)
      {
        weta += E77[i][j] * di * di;
        wphi += E77[i][j] * dj * dj;
        weta_cog += E77[i][j] * di_float * di_float;
        wphi_cog += E77[i][j] * dj_float * dj_float;
        Eetaphi += E77[i][j];
        if (i != 3 || j != 3)
        {
          weta_cogx += E77[i][j] * di_float * di_float;
          wphi_cogx += E77[i][j] * dj_float * dj_float;
        }
      }

      e77 += E77[i][j];
      if (di <= 1 && (dj == 0 || j == (3 + signphi)))
      {
        w32 += E77[i][j] * (i - 3) * (i - 3);
        e32 += E77[i][j];
      }
      if (di <= 2 && (dj == 0 || j == (3 + signphi)))
      {
        w52 += E77[i][j] * (i - 3) * (i - 3);
        e52 += E77[i][j];
      }
      if (di <= 3 && (dj == 0 || j == (3 + signphi)))
      {
        w72 += E77[i][j] * (i - 3) * (i - 3);
        e72 += E77[i][j];
      }

      if (di <= 0 && dj <= 1)
      {
        e13 += E77[i][j];
      }
      if (di <= 0 && dj <= 2)
      {
        e15 += E77[i][j];
      }
      if (di <= 0 && dj <= 3)
      {
        e17 += E77[i][j];
      }
      if (di <= 1 && dj <= 0)
      {
        e31 += E77[i][j];
      }
      if (di <= 2 && dj <= 0)
      {
        e51 += E77[i][j];
      }
      if (di <= 3 && dj <= 0)
      {
        e71 += E77[i][j];
      }
      if (di <= 1 && dj <= 1)
      {
        e33 += E77[i][j];
      }
      if (di <= 1 && dj <= 2)
      {
        e35 += E77[i][j];
      }
      if (di <= 1 && dj <= 3)
      {
        e37 += E77[i][j];
      }
      if (di <= 2 && dj <= 1)
      {
        e53 += E77[i][j];
      }
      if (di <= 3 && dj <= 1)
      {
        e73 += E77[i][j];
      }
      if (di <= 2 && dj <= 2)
      {
        e55 += E77[i][j];
      }
      if (di <= 2 && dj <= 3)
      {
        e57 += E77[i][j];
      }
      if (di <= 3 && dj <= 2)
      {
        e75 += E77[i][j];
      }
    }
  }

  if (Eetaphi > 0)
  {
    weta /= Eetaphi;
    wphi /= Eetaphi;
    weta_cog /= Eetaphi;
    wphi_cog /= Eetaphi;
    weta_cogx /= Eetaphi;
    wphi_cogx /= Eetaphi;
  }
  if (e32 > 0)
  {
    w32 /= e32;
  }
  if (e52 > 0)
  {
    w52 /= e52;
  }
  if (e72 > 0)
  {
    w72 /= e72;
  }

  photon->set_shower_shape_parameter("et1", showershape[0]);
  photon->set_shower_shape_parameter("et2", showershape[1]);
  photon->set_shower_shape_parameter("et3", showershape[2]);
  photon->set_shower_shape_parameter("et4", showershape[3]);
    photon->set_shower_shape_parameter("e11", e11);

    float e22 = -999.0F;
    if (showershape.size() >= 12)
    {
      e22 = showershape[8] + showershape[9] + showershape[10] + showershape[11];
    }
    photon->set_shower_shape_parameter("e22", e22);

    photon->set_shower_shape_parameter("e33", e33);
  photon->set_shower_shape_parameter("e55", e55);
  photon->set_shower_shape_parameter("e77", e77);
  photon->set_shower_shape_parameter("e13", e13);
  photon->set_shower_shape_parameter("e15", e15);
  photon->set_shower_shape_parameter("e17", e17);
  photon->set_shower_shape_parameter("e31", e31);
  photon->set_shower_shape_parameter("e51", e51);
  photon->set_shower_shape_parameter("e71", e71);
  photon->set_shower_shape_parameter("e35", e35);
  photon->set_shower_shape_parameter("e37", e37);
  photon->set_shower_shape_parameter("e53", e53);
  photon->set_shower_shape_parameter("e73", e73);
  photon->set_shower_shape_parameter("e57", e57);
  photon->set_shower_shape_parameter("e75", e75);
  photon->set_shower_shape_parameter("weta", weta);
  photon->set_shower_shape_parameter("wphi", wphi);
  photon->set_shower_shape_parameter("weta_cog", weta_cog);
  photon->set_shower_shape_parameter("wphi_cog", wphi_cog);
  photon->set_shower_shape_parameter("weta_cogx", weta_cogx);
  photon->set_shower_shape_parameter("wphi_cogx", wphi_cogx);
  photon->set_shower_shape_parameter("detamax", detamax);
  photon->set_shower_shape_parameter("dphimax", dphimax);
  photon->set_shower_shape_parameter("nsaturated", nsaturated);
  photon->set_shower_shape_parameter("e32", e32);
  photon->set_shower_shape_parameter("e52", e52);
  photon->set_shower_shape_parameter("e72", e72);
  photon->set_shower_shape_parameter("w32", w32);
  photon->set_shower_shape_parameter("w52", w52);
  photon->set_shower_shape_parameter("w72", w72);
  photon->set_shower_shape_parameter("cluster_eta", cluster_eta);
  photon->set_shower_shape_parameter("cluster_phi", cluster_phi);
  photon->set_shower_shape_parameter("mean_time", clusteravgtime);
  photon->set_shower_shape_parameter("detacog", detacog);
  photon->set_shower_shape_parameter("dphicog", dphicog);
  photon->set_shower_shape_parameter("drad", drad);

  // HCAL info
  std::vector<int> ihcal_tower = find_closest_hcal_tower(cluster_eta, cluster_phi, m_geomIH, m_ihcal_tower_container, 0.0, true);
  std::vector<int> ohcal_tower = find_closest_hcal_tower(cluster_eta, cluster_phi, m_geomOH, m_ohcal_tower_container, 0.0, false);

  float ihcal_et = 0;
  float ohcal_et = 0;
  float ihcal_et22 = 0;
  float ohcal_et22 = 0;
  float ihcal_et33 = 0;
  float ohcal_et33 = 0;

  int ihcal_ieta = ihcal_tower[0];
  int ihcal_iphi = ihcal_tower[1];
  float ihcalEt33[3][3] = {{0.0F}};

  int ohcal_ieta = ohcal_tower[0];
  int ohcal_iphi = ohcal_tower[1];
  float ohcalEt33[3][3] = {{0.0F}};

  for (int ieta_h = ihcal_ieta - 1; ieta_h <= ihcal_ieta + 1; ieta_h++)
  {
    for (int iphi_h = ihcal_iphi - 1; iphi_h <= ihcal_iphi + 1; iphi_h++)
    {
      int temp_ieta = ieta_h;
      int temp_iphi = iphi_h;
      shift_tower_index(temp_ieta, temp_iphi, 24, 64);
      if (temp_ieta < 0)
      {
        continue;
      }

      unsigned int towerinfokey = TowerInfoDefs::encode_hcal(temp_ieta, temp_iphi);
      TowerInfo* towerinfo = m_ihcal_tower_container->get_tower_at_key(towerinfokey);
      if (towerinfo && towerinfo->get_isGood())
      {
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, temp_ieta, temp_iphi);
        RawTowerGeom* tower_geom = m_geomIH->get_tower_geometry(key);
        if (tower_geom)
        {
          float energy = towerinfo->get_energy();
          float eta = getTowerEta(tower_geom, 0, 0, 0);
          float sintheta = 1.0 / std::cosh(eta);
          float Et = energy * sintheta;
          ihcalEt33[ieta_h - ihcal_ieta + 1][iphi_h - ihcal_iphi + 1] = Et;
        }
      }
    }
  }

  for (int ieta_h = ohcal_ieta - 1; ieta_h <= ohcal_ieta + 1; ieta_h++)
  {
    for (int iphi_h = ohcal_iphi - 1; iphi_h <= ohcal_iphi + 1; iphi_h++)
    {
      int temp_ieta = ieta_h;
      int temp_iphi = iphi_h;
      shift_tower_index(temp_ieta, temp_iphi, 24, 64);
      if (temp_ieta < 0)
      {
        continue;
      }

      unsigned int towerinfokey = TowerInfoDefs::encode_hcal(temp_ieta, temp_iphi);
      TowerInfo* towerinfo = m_ohcal_tower_container->get_tower_at_key(towerinfokey);
      if (towerinfo && towerinfo->get_isGood())
      {
        const RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, temp_ieta, temp_iphi);
        RawTowerGeom* tower_geom = m_geomOH->get_tower_geometry(key);
        if (tower_geom)
        {
          float energy = towerinfo->get_energy();
          float eta = getTowerEta(tower_geom, 0, 0, 0);
          float sintheta = 1.0 / std::cosh(eta);
          float Et = energy * sintheta;
          ohcalEt33[ieta_h - ohcal_ieta + 1][iphi_h - ohcal_iphi + 1] = Et;
        }
      }
    }
  }

  ihcal_et = ihcalEt33[1][1];
  ohcal_et = ohcalEt33[1][1];

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      ihcal_et33 += ihcalEt33[i][j];
      ohcal_et33 += ohcalEt33[i][j];
      if (i == 1 || j == 1 + ihcal_tower[2])
      {
        if (j == 1 || i == 1 + ihcal_tower[3])
        {
          ihcal_et22 += ihcalEt33[i][j];
        }
      }
      if (i == 1 || j == 1 + ohcal_tower[2])
      {
        if (j == 1 || i == 1 + ohcal_tower[3])
        {
          ohcal_et22 += ohcalEt33[i][j];
        }
      }
    }
  }

  photon->set_shower_shape_parameter("ihcal_et", ihcal_et);
  photon->set_shower_shape_parameter("ohcal_et", ohcal_et);
  photon->set_shower_shape_parameter("ihcal_et22", ihcal_et22);
  photon->set_shower_shape_parameter("ohcal_et22", ohcal_et22);
  photon->set_shower_shape_parameter("ihcal_et33", ihcal_et33);
  photon->set_shower_shape_parameter("ohcal_et33", ohcal_et33);
  photon->set_shower_shape_parameter("ihcal_ieta", ihcal_ieta);
  photon->set_shower_shape_parameter("ihcal_iphi", ihcal_iphi);
  photon->set_shower_shape_parameter("ohcal_ieta", ohcal_ieta);
  photon->set_shower_shape_parameter("ohcal_iphi", ohcal_iphi);

  float E = photon->get_energy();
  float ET = E / std::cosh(cluster_eta);

  auto compute_layer_iso = [&](RawTowerDefs::CalorimeterId calo_id, float radius)
  {
    TowerInfoContainer* container = nullptr;
    RawTowerGeomContainer* geom = nullptr;
    if (calo_id == RawTowerDefs::CalorimeterId::CEMC)
    {
      container = m_emc_tower_container;
      geom = m_geomEM;
    }
    else if (calo_id == RawTowerDefs::CalorimeterId::HCALIN)
    {
      container = m_ihcal_tower_container;
      geom = m_geomIH;
    }
    else
    {
      container = m_ohcal_tower_container;
      geom = m_geomOH;
    }
    return calculate_layer_et(cluster_eta, cluster_phi, radius, container, geom, calo_id, m_vertex);
  };

  const float emcal_et_04 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.4);
  const float ihcal_et_04 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALIN, 0.4);
  const float ohcal_et_04 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALOUT, 0.4);

  const float emcal_et_03 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.3);
  const float ihcal_et_03 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALIN, 0.3);
  const float ohcal_et_03 = compute_layer_iso(RawTowerDefs::CalorimeterId::HCALOUT, 0.3);

  const float emcal_et_02 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.2);
  const float emcal_et_01 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.1);
  const float emcal_et_005 = compute_layer_iso(RawTowerDefs::CalorimeterId::CEMC, 0.05);

    // -----------------------------
    // Build iso values (same as before)
    // -----------------------------
    const float iso_04_emcal  = emcal_et_04  - ET;
    const float iso_03_emcal  = emcal_et_03  - ET;
    const float iso_02_emcal  = emcal_et_02  - ET;
    const float iso_01_emcal  = emcal_et_01  - ET;
    const float iso_005_emcal = emcal_et_005 - ET;

    // Total iso (what RecoilJets effectively uses when it adds EMCal + HCal pieces)
    const float iso03_total = iso_03_emcal + ihcal_et_03 + ohcal_et_03;
    const float iso04_total = iso_04_emcal + ihcal_et_04 + ohcal_et_04;

    // -----------------------------
    // NEW: print ONLY when negative iso appears
    // -----------------------------
    if (Verbosity() >= 2 && (iso03_total < -0.05f || iso04_total < -0.05f))
    {
      static int s_negPrint = 0;

      std::cout << Name()
                << ": NEG ISO detected"
                << " | vertex_z=" << m_vertex
                << " | eta=" << cluster_eta
                << " phi=" << cluster_phi
                << " | E(rc)=" << rc->get_energy()
                << " ET=" << ET
                << " | emcal_et_03=" << emcal_et_03
                << " iso_03_emcal=" << iso_03_emcal
                << " (ih=" << ihcal_et_03 << ", oh=" << ohcal_et_03 << ", total=" << iso03_total << ")"
                << " | emcal_et_04=" << emcal_et_04
                << " iso_04_emcal=" << iso_04_emcal
                << " (ih=" << ihcal_et_04 << ", oh=" << ohcal_et_04 << ", total=" << iso04_total << ")"
                << std::endl;

      // Extra, cheap sanity: do the *cluster's own towers* have TowerInfo energy consistent with rc->get_energy()?
      // (This catches the "TowerInfo energy is 0 / mismatched container" problem immediately.)
        {
          // ================================================================
          // (A) Cluster-tower sanity: compare RawCluster energy vs TowerInfo energies
          // ================================================================
          const RawCluster::TowerMap& tower_map_dbg = rc->get_towermap();

          int   nTI_found   = 0;
          int   nTI_nonzero = 0;
          float sumE_TI     = 0.0f;
          float sumEt_TI    = 0.0f;

          for (const auto& it : tower_map_dbg)
          {
            RawTowerDefs::keytype tower_key = it.first;
            int ieta = RawTowerDefs::decode_index1(tower_key);
            int iphi = RawTowerDefs::decode_index2(tower_key);

            unsigned int towerinfokey = TowerInfoDefs::encode_emcal(ieta, iphi);
            TowerInfo* ti = m_emc_tower_container ? m_emc_tower_container->get_tower_at_key(towerinfokey) : nullptr;
            if (!ti) continue;

            ++nTI_found;

            const float e = ti->get_energy();
            if (e > 0.0f) ++nTI_nonzero;
            sumE_TI += e;

            RawTowerDefs::keytype geom_key =
              RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
            RawTowerGeom* tg = m_geomEM ? m_geomEM->get_tower_geometry(geom_key) : nullptr;
            if (tg)
            {
              const double tEta = getTowerEta(tg, 0, 0, m_vertex);
              if (std::isfinite(tEta)) sumEt_TI += (e / std::cosh(tEta));
            }
          }

          std::cout << "  [cluster TI sanity]"
                    << " towers_in_cluster=" << tower_map_dbg.size()
                    << " TI_found=" << nTI_found
                    << " TI_nonzero=" << nTI_nonzero
                    << " sumE_TI=" << sumE_TI
                    << " sumEt_TI=" << sumEt_TI
                    << " rc_E=" << rc->get_energy()
                    << std::endl;

          // ================================================================
          // (B) Cone-sum audit: reproduce calculate_layer_et selection and count DROP reasons
          //     Only meaningful for EMCal here because NEG iso is dominated by (emcal_et - ET).
          // ================================================================
          auto auditCone = [&](float Rcone)
          {
            if (!m_emc_tower_container || !m_geomEM)
            {
              std::cout << "  [cone audit] R=" << Rcone << " : missing TOWERINFO_CALIB_CEMC or TOWERGEOM_CEMC\n";
              return;
            }

            const unsigned int ntowers = m_emc_tower_container->size();

            // counters
            unsigned int nNull      = 0;
            unsigned int nNotGood   = 0;
            unsigned int nGeomMiss  = 0;
            unsigned int nOutR      = 0;
            unsigned int nBelowThr  = 0;
            unsigned int nKept      = 0;

            double sumEtKept = 0.0;

            struct Row
            {
              int ieta = -1;
              int iphi = -1;
              double dR = 0.0;
              double e  = 0.0;
              double et = 0.0;
              int isGood = 0;
              const char* why = "";
            };

            std::vector<Row> droppedInCone;
            std::vector<Row> keptInCone;
            droppedInCone.reserve(64);
            keptInCone.reserve(64);

            for (unsigned int ch = 0; ch < ntowers; ++ch)
            {
              TowerInfo* tower = m_emc_tower_container->get_tower_at_channel(ch);
              if (!tower) { ++nNull; continue; }

              if (!tower->get_isGood())
              {
                // still need dR to decide if it's "relevant" (in cone)
                // but dR requires geometry; if geometry missing, count there.
                // We handle dR classification below once geometry exists.
              }

              const unsigned int tkey = m_emc_tower_container->encode_key(ch);
              const int ieta = m_emc_tower_container->getTowerEtaBin(tkey);
              const int iphi = m_emc_tower_container->getTowerPhiBin(tkey);

              RawTowerDefs::keytype geom_key =
                RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
              RawTowerGeom* tg = m_geomEM->get_tower_geometry(geom_key);
              if (!tg) { ++nGeomMiss; continue; }

              const double tEta = getTowerEta(tg, 0, 0, m_vertex);
              const double tPhi = tg->get_phi();
              if (!std::isfinite(tEta) || !std::isfinite(tPhi)) continue;

              const double dRval = deltaR(cluster_eta, cluster_phi, tEta, tPhi);

              // outside cone
              if (dRval >= Rcone) { ++nOutR; continue; }

              // now inside cone: classify pass/fail reasons exactly like calculate_layer_et
              const float e = tower->get_energy();
              const int isGood = tower->get_isGood() ? 1 : 0;

              if (!tower->get_isGood())
              {
                ++nNotGood;
                droppedInCone.push_back({ieta, iphi, dRval, e, 0.0, isGood, "!isGood"});
                continue;
              }

              if (e <= m_shape_min_tower_E)
              {
                ++nBelowThr;
                droppedInCone.push_back({ieta, iphi, dRval, e, 0.0, isGood, "E<=min"});
                continue;
              }

              const double et = e / std::cosh(tEta);
              ++nKept;
              sumEtKept += et;
              keptInCone.push_back({ieta, iphi, dRval, e, et, isGood, "KEPT"});
            }

            // Sort kept towers by Et descending; dropped by E descending (just to show informative ones)
            std::sort(keptInCone.begin(), keptInCone.end(),
                      [](const Row& a, const Row& b){ return a.et > b.et; });
            std::sort(droppedInCone.begin(), droppedInCone.end(),
                      [](const Row& a, const Row& b){ return a.e > b.e; });

            std::cout << "\n  [cone audit] EMCal R=" << std::fixed << std::setprecision(2) << Rcone
                      << "  seed(eta,phi)=(" << std::setprecision(3) << cluster_eta
                      << "," << std::setprecision(3) << cluster_phi << ")"
                      << "  vertex_z=" << std::setprecision(2) << m_vertex
                      << "\n"
                      << "    ntowers=" << ntowers
                      << " | null=" << nNull
                      << " | geomMissing=" << nGeomMiss
                      << " | outR=" << nOutR
                      << " | inR_notGood=" << nNotGood
                      << " | inR_E<=min(" << m_shape_min_tower_E << ")=" << nBelowThr
                      << " | inR_KEPT=" << nKept
                      << " | sumEtKept=" << std::setprecision(3) << sumEtKept
                      << "\n";

            auto printTable = [&](const char* title, const std::vector<Row>& rows, std::size_t maxRows)
            {
              std::cout << "    " << title << " (showing up to " << maxRows << ")\n";
              std::cout << "    "
                        << std::setw(6)  << "ieta"
                        << std::setw(6)  << "iphi"
                        << std::setw(10) << "dR"
                        << std::setw(12) << "E"
                        << std::setw(12) << "Et"
                        << std::setw(8)  << "isGood"
                        << "  why\n";
              std::cout << "    ---------------------------------------------------------------\n";

              const std::size_t n = std::min<std::size_t>(rows.size(), maxRows);
              for (std::size_t i = 0; i < n; ++i)
              {
                const Row& r = rows[i];
                std::cout << "    "
                          << std::setw(6)  << r.ieta
                          << std::setw(6)  << r.iphi
                          << std::setw(10) << std::fixed << std::setprecision(4) << r.dR
                          << std::setw(12) << std::fixed << std::setprecision(4) << r.e
                          << std::setw(12) << std::fixed << std::setprecision(4) << r.et
                          << std::setw(8)  << r.isGood
                          << "  " << r.why
                          << "\n";
              }
            };

            printTable("KEPT towers inside cone (highest Et first)", keptInCone, 12);
            if (!droppedInCone.empty())
              printTable("DROPPED towers inside cone (highest E first)", droppedInCone, 12);
          };

          // Audit the two radii you actually use for iso fields
          auditCone(0.30f);
          auditCone(0.40f);
        }

      // Print only first ~15 occurrences to avoid log spam
      ++s_negPrint;
      if (s_negPrint >= 15)
      {
        std::cout << Name()
                  << ": NEG ISO print limit reached (15). "
                  << "Raise/adjust this cap if you need more examples."
                  << std::endl;
      }
    }

    // -----------------------------
    // Store to PhotonClusterv1 (same names as before)
    // -----------------------------
    photon->set_shower_shape_parameter("iso_04_emcal",  iso_04_emcal);
    photon->set_shower_shape_parameter("iso_04_hcalin", ihcal_et_04);
    photon->set_shower_shape_parameter("iso_04_hcalout",ohcal_et_04);

    photon->set_shower_shape_parameter("iso_03_emcal",  iso_03_emcal);
    photon->set_shower_shape_parameter("iso_03_hcalin", ihcal_et_03);
    photon->set_shower_shape_parameter("iso_03_hcalout",ohcal_et_03);

    photon->set_shower_shape_parameter("iso_02_emcal",  iso_02_emcal);
    photon->set_shower_shape_parameter("iso_01_emcal",  iso_01_emcal);
    photon->set_shower_shape_parameter("iso_005_emcal", iso_005_emcal);
}

double PhotonClusterBuilder::getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz)
{
  if (!tower_geom)
  {
    return -9999;
  }
  if (vx == 0 && vy == 0 && vz == 0)
  {
    return tower_geom->get_eta();
  }

  double radius = sqrt((tower_geom->get_center_x() - vx) * (tower_geom->get_center_x() - vx) + (tower_geom->get_center_y() - vy) * (tower_geom->get_center_y() - vy));
  double theta = atan2(radius, tower_geom->get_center_z() - vz);
  return -log(tan(theta / 2.));
}

std::vector<int> PhotonClusterBuilder::find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer* geom, TowerInfoContainer* towerContainer, float vertex_z, bool isihcal)
{
  int matchedieta = -1;
  int matchediphi = -1;
  double matchedeta = -999;
  double matchedphi = -999;

  if (!geom || !towerContainer)
  {
    return {-1, -1, 0, 0};
  }

  unsigned int ntowers = towerContainer->size();
  float minR = 999;

  for (unsigned int channel = 0; channel < ntowers; channel++)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower)
    {
      continue;
    }

    unsigned int towerkey = towerContainer->encode_key(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);

    RawTowerDefs::keytype key = RawTowerDefs::encode_towerid(isihcal ? RawTowerDefs::CalorimeterId::HCALIN : RawTowerDefs::CalorimeterId::HCALOUT, ieta, iphi);
    RawTowerGeom* tower_geom = geom->get_tower_geometry(key);
    if (!tower_geom)
    {
      continue;
    }

    double this_phi = tower_geom->get_phi();
    double this_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double dR_val = deltaR(eta, phi, this_eta, this_phi);
    if (dR_val < minR)
    {
      minR = dR_val;
      matchedieta = ieta;
      matchediphi = iphi;
      matchedeta = this_eta;
      matchedphi = this_phi;
    }
  }

  float deta = eta - matchedeta;
  float dphi_val = phi - matchedphi;
  if (dphi_val > M_PI)
  {
    dphi_val -= 2 * M_PI;
  }
  if (dphi_val < -M_PI)
  {
    dphi_val += 2 * M_PI;
  }

  int dphisign = (dphi_val > 0) ? 1 : -1;
  int detasign = (deta > 0) ? 1 : -1;

  return {matchedieta, matchediphi, detasign, dphisign};
}

float PhotonClusterBuilder::calculate_layer_et(float seed_eta, float seed_phi, float radius, TowerInfoContainer* towerContainer, RawTowerGeomContainer* geomContainer, RawTowerDefs::CalorimeterId calo_id, float vertex_z)
{
  if (!towerContainer || !geomContainer)
  {
    return std::numeric_limits<float>::quiet_NaN();
  }

  float layer_et = 0.0;
  const unsigned int ntowers = towerContainer->size();
  for (unsigned int channel = 0; channel < ntowers; ++channel)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower || !tower->get_isGood())
    {
      continue;
    }

    unsigned int towerkey = towerContainer->encode_key(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);

    RawTowerDefs::keytype geom_key = RawTowerDefs::encode_towerid(calo_id, ieta, iphi);
    RawTowerGeom* tower_geom = geomContainer->get_tower_geometry(geom_key);
    if (!tower_geom)
    {
      continue;
    }

    double tower_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
    double tower_phi = tower_geom->get_phi();

    if (deltaR(seed_eta, seed_phi, tower_eta, tower_phi) >= radius)
    {
      continue;
    }

    float energy = tower->get_energy();
    if (energy <= m_shape_min_tower_E)
    {
      continue;
    }

    float et = energy / std::cosh(tower_eta);
    layer_et += et;
  }

  return layer_et;
}

double PhotonClusterBuilder::deltaR(double eta1, double phi1, double eta2, double phi2)
{
  double dphi = phi1 - phi2;
  while (dphi > M_PI)
  {
    dphi -= 2 * M_PI;
  }
  while (dphi <= -M_PI)
  {
    dphi += 2 * M_PI;
  }
  return sqrt(pow(eta1 - eta2, 2) + pow(dphi, 2));
}
