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
#include <cstdlib>
#include <iostream>
#include <iomanip>   // NEW: for std::setw / std::setprecision in debug tables
#include <set>
#include <sstream>
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

    m_emc_tower_container = findNode::getClass<TowerInfoContainer>(topNode, m_emc_tower_node);
    if (!m_emc_tower_container)
    {
      std::cerr << Name() << ": could not find TowerInfoContainer node '" << m_emc_tower_node << "'" << std::endl;
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

    m_ihcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, m_ihcal_tower_node);
    if (!m_ihcal_tower_container)
    {
      std::cerr << Name() << ": could not find TowerInfoContainer node '" << m_ihcal_tower_node << "'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_geomIH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if (!m_geomIH)
    {
      std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALIN'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_ohcal_tower_container = findNode::getClass<TowerInfoContainer>(topNode, m_ohcal_tower_node);
    if (!m_ohcal_tower_container)
    {
      std::cerr << Name() << ": could not find TowerInfoContainer node '" << m_ohcal_tower_node << "'" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_geomOH = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if (!m_geomOH)
      {
        std::cerr << Name() << ": could not find RawTowerGeomContainer node 'TOWERGEOM_HCALOUT'" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    // ------------------------------------------------------------
    // Optional: Au+Au UE-subtracted isolation (cone ET sums)
    //   - Keep m_emc_tower_container / m_geomEM for shower-shapes (96x256 CEMC)
    //   - For isolation only, switch to *_SUB1 nodes when requested
    //   - NOTE: CEMC_RETOWER_SUB1 is indexed like HCALIN (24x64), so geometry/id must use HCALIN
    // ------------------------------------------------------------
    if (m_is_auau)
    {
      const std::string cemcIsoNode  = m_tower_node_prefix + "_CEMC_RETOWER_SUB1";
      const std::string ihcalIsoNode = m_tower_node_prefix + "_HCALIN_SUB1";
      const std::string ohcalIsoNode = m_tower_node_prefix + "_HCALOUT_SUB1";

      m_emc_tower_container_iso = findNode::getClass<TowerInfoContainer>(topNode, cemcIsoNode);
      m_ihcal_tower_container_iso = findNode::getClass<TowerInfoContainer>(topNode, ihcalIsoNode);
      m_ohcal_tower_container_iso = findNode::getClass<TowerInfoContainer>(topNode, ohcalIsoNode);

      // Retowered CEMC uses HCALIN geometry
      m_geomEM_iso = m_geomIH;

      if (!m_emc_tower_container_iso || !m_ihcal_tower_container_iso || !m_ohcal_tower_container_iso || !m_geomEM_iso)
      {
        std::cerr << Name()
                  << ": AuAu UE-subtracted isolation requested but required nodes are missing.\n"
                  << "  expected: " << cemcIsoNode << ", " << ihcalIsoNode << ", " << ohcalIsoNode
                  << " (and TOWERGEOM_HCALIN)\n";
        return Fun4AllReturnCodes::ABORTRUN;
      }

      if (Verbosity() > 0)
      {
        std::cout << Name() << ": AuAu UE-subtracted iso enabled"
                  << " | prefix=" << m_tower_node_prefix
                  << " | cemc_iso=" << cemcIsoNode
                  << " | ihcal_iso=" << ihcalIsoNode
                  << " | ohcal_iso=" << ohcalIsoNode
                  << std::endl;
      }
    }
    else
    {
      m_emc_tower_container_iso = nullptr;
      m_ihcal_tower_container_iso = nullptr;
      m_ohcal_tower_container_iso = nullptr;
      m_geomEM_iso = nullptr;
    }

    m_iso_audit_mode = false;
    if (const char* env = std::getenv("RJ_ISO_AUDIT_MODE"))
    {
      m_iso_audit_mode = (std::atoi(env) != 0);
    }

    m_iso_audit_summary_every_events = 1000;
    if (const char* env_builder = std::getenv("RJ_ISO_AUDIT_BUILDER_SUMMARY_EVERY_EVENTS"))
    {
      const int v = std::atoi(env_builder);
      if (v > 0) m_iso_audit_summary_every_events = v;
    }
    else if (const char* env_skip = std::getenv("RJ_ISO_AUDIT_SKIP_SUMMARY_EVERY_EVENTS"))
    {
      const int v = std::atoi(env_skip);
      if (v > 0) m_iso_audit_summary_every_events = v;
    }

    m_evt_seen = 0;
    m_evt_missing_rawclusters = 0;
    m_evt_vertex_from_mbd = 0;
    m_evt_vertex_from_global = 0;
    m_evt_no_finite_vertex = 0;
    m_evt_skip_vz = 0;
    m_evt_zero_input_clusters = 0;
    m_evt_zero_pass_et_clusters = 0;
    m_evt_zero_built_after_pass_et = 0;
    m_evt_built_photons = 0;
    m_clusters_seen_total = 0;
    m_clusters_pass_et_total = 0;
    m_photons_built_total = 0;
    m_mean_time_default_total = 0;
    m_nonfinite_eta_total = 0;
    m_nonfinite_phi_total = 0;
    m_nonfinite_et_total = 0;
    m_audit_last_summary = AuditSnapshot{};

    if (m_iso_audit_mode)
    {
      std::cout << Name()
                << ": [InitRun][audit]"
                << " input=" << m_input_cluster_node
                << " output=" << m_output_photon_node
                << " ETthr=" << m_min_cluster_et
                << " shapeTowerMinE=" << m_shape_min_tower_E
                << " isoTowerPolicy=no_one_sided_cut"
                << " useVzCut=" << (m_use_vz_cut ? "true" : "false")
                << " vzCut=" << m_vz_cut_cm
                << " isAuAu=" << (m_is_auau ? "true" : "false")
                << std::endl;

      std::cout << Name()
                << ": [InitRun][audit]"
                << " towers(EM/IH/OH)="
                << m_emc_tower_node << " / "
                << m_ihcal_tower_node << " / "
                << m_ohcal_tower_node
                << " | towerPrefix=" << m_tower_node_prefix
                << " | cadence=" << m_iso_audit_summary_every_events
                << std::endl;
    }

    CreateNodes(topNode);
    m_audit_last_summary = make_audit_snapshot();
    return Fun4AllReturnCodes::EVENT_OK;
}

PhotonClusterBuilder::AuditSnapshot PhotonClusterBuilder::make_audit_snapshot() const
{
  AuditSnapshot s{};
  s.evt_seen = m_evt_seen;
  s.evt_missing_rawclusters = m_evt_missing_rawclusters;
  s.evt_vertex_from_mbd = m_evt_vertex_from_mbd;
  s.evt_vertex_from_global = m_evt_vertex_from_global;
  s.evt_no_finite_vertex = m_evt_no_finite_vertex;
  s.evt_skip_vz = m_evt_skip_vz;
  s.evt_zero_input_clusters = m_evt_zero_input_clusters;
  s.evt_zero_pass_et_clusters = m_evt_zero_pass_et_clusters;
  s.evt_zero_built_after_pass_et = m_evt_zero_built_after_pass_et;
  s.evt_built_photons = m_evt_built_photons;
  s.clusters_seen_total = m_clusters_seen_total;
  s.clusters_pass_et_total = m_clusters_pass_et_total;
  s.photons_built_total = m_photons_built_total;
  s.mean_time_default_total = m_mean_time_default_total;
  s.nonfinite_eta_total = m_nonfinite_eta_total;
  s.nonfinite_phi_total = m_nonfinite_phi_total;
  s.nonfinite_et_total = m_nonfinite_et_total;
  return s;
}

void PhotonClusterBuilder::print_audit_summary(bool force)
{
  if (!m_iso_audit_mode) return;

  const unsigned long long cadence =
      (m_iso_audit_summary_every_events > 0 ? static_cast<unsigned long long>(m_iso_audit_summary_every_events) : 1000ULL);

  if (!force)
  {
    if (m_evt_seen == 0ULL) return;
    if ((m_evt_seen % cadence) != 0ULL) return;
  }

  const AuditSnapshot cur = make_audit_snapshot();
  const AuditSnapshot prev = m_audit_last_summary;

  auto diff = [](unsigned long long a, unsigned long long b) -> unsigned long long
  {
    return (a >= b ? (a - b) : 0ULL);
  };

  const unsigned long long win_evt_seen = diff(cur.evt_seen, prev.evt_seen);
  const unsigned long long win_evt_missing_rawclusters = diff(cur.evt_missing_rawclusters, prev.evt_missing_rawclusters);
  const unsigned long long win_evt_vertex_from_mbd = diff(cur.evt_vertex_from_mbd, prev.evt_vertex_from_mbd);
  const unsigned long long win_evt_vertex_from_global = diff(cur.evt_vertex_from_global, prev.evt_vertex_from_global);
  const unsigned long long win_evt_no_finite_vertex = diff(cur.evt_no_finite_vertex, prev.evt_no_finite_vertex);
  const unsigned long long win_evt_skip_vz = diff(cur.evt_skip_vz, prev.evt_skip_vz);
  const unsigned long long win_evt_zero_input_clusters = diff(cur.evt_zero_input_clusters, prev.evt_zero_input_clusters);
  const unsigned long long win_evt_zero_pass_et_clusters = diff(cur.evt_zero_pass_et_clusters, prev.evt_zero_pass_et_clusters);
  const unsigned long long win_evt_zero_built_after_pass_et = diff(cur.evt_zero_built_after_pass_et, prev.evt_zero_built_after_pass_et);
  const unsigned long long win_evt_built_photons = diff(cur.evt_built_photons, prev.evt_built_photons);

  const unsigned long long win_clusters_seen_total = diff(cur.clusters_seen_total, prev.clusters_seen_total);
  const unsigned long long win_clusters_pass_et_total = diff(cur.clusters_pass_et_total, prev.clusters_pass_et_total);
  const unsigned long long win_photons_built_total = diff(cur.photons_built_total, prev.photons_built_total);
  const unsigned long long win_mean_time_default_total = diff(cur.mean_time_default_total, prev.mean_time_default_total);
  const unsigned long long win_nonfinite_eta_total = diff(cur.nonfinite_eta_total, prev.nonfinite_eta_total);
  const unsigned long long win_nonfinite_phi_total = diff(cur.nonfinite_phi_total, prev.nonfinite_phi_total);
  const unsigned long long win_nonfinite_et_total = diff(cur.nonfinite_et_total, prev.nonfinite_et_total);

  const char* dominant_label = "none";
  unsigned long long dominant_count = 0;

  auto consider = [&](unsigned long long count, const char* label)
  {
    if (count > dominant_count)
    {
      dominant_count = count;
      dominant_label = label;
    }
  };

  consider(win_evt_no_finite_vertex, "noFiniteRecoVertex");
  consider(win_evt_skip_vz, "skipVz");
  consider(win_evt_zero_input_clusters, "zeroInputClusters");
  consider(win_evt_zero_pass_et_clusters, "zeroPassET");
  consider(win_evt_zero_built_after_pass_et, "zeroBuiltAfterPassET");

  std::cout << Name()
            << ": [" << (force ? "FINAL" : "SUMMARY") << "]"
            << " evt=" << m_evt_seen
            << " windowEvents=" << win_evt_seen
            << " cadence=" << m_iso_audit_summary_every_events
            << " | input=" << m_input_cluster_node
            << " | output=" << m_output_photon_node
            << " | ETthr=" << m_min_cluster_et
            << " | vzCut=" << m_vz_cut_cm
            << std::endl;

  std::cout << "  Δevents:"
            << " rawclustersMissing=" << win_evt_missing_rawclusters
            << " vertex(MBD/Global)=" << win_evt_vertex_from_mbd << "/" << win_evt_vertex_from_global
            << " noFiniteRecoVertex=" << win_evt_no_finite_vertex
            << " skipVz=" << win_evt_skip_vz
            << " zeroInputClusters=" << win_evt_zero_input_clusters
            << " zeroPassET=" << win_evt_zero_pass_et_clusters
            << " zeroBuiltAfterPassET=" << win_evt_zero_built_after_pass_et
            << " built>0=" << win_evt_built_photons
            << std::endl;

  std::cout << "  Δclusters:"
            << " seen=" << win_clusters_seen_total
            << " passET=" << win_clusters_pass_et_total
            << " built=" << win_photons_built_total
            << " nonFinite(eta,phi,ET)=("
            << win_nonfinite_eta_total << ","
            << win_nonfinite_phi_total << ","
            << win_nonfinite_et_total << ")"
            << " mean_time=-999=" << win_mean_time_default_total
            << std::endl;

  if (dominant_count > 0ULL)
  {
    std::cout << "  dominant zero-photon cause in window: "
              << dominant_label << "=" << dominant_count;
    if (std::string(dominant_label) == "zeroPassET")
    {
      std::cout << " (all input clusters stayed below the ET threshold)";
    }
    std::cout << std::endl;
  }
  else
  {
    std::cout << "  dominant zero-photon cause in window: none" << std::endl;
  }

  std::cout << "  cumulative:"
            << " events=" << cur.evt_seen
            << " rawclustersMissing=" << cur.evt_missing_rawclusters
            << " skipVz=" << cur.evt_skip_vz
            << " zeroInputClusters=" << cur.evt_zero_input_clusters
            << " zeroPassET=" << cur.evt_zero_pass_et_clusters
            << " zeroBuiltAfterPassET=" << cur.evt_zero_built_after_pass_et
            << " built>0=" << cur.evt_built_photons
            << " clusters=" << cur.clusters_seen_total
            << " passET=" << cur.clusters_pass_et_total
            << " built=" << cur.photons_built_total
            << std::endl;

  if (force && win_evt_seen == 0ULL)
  {
    std::cout << "  note: no new events were processed since the last periodic builder summary." << std::endl;
  }

  m_audit_last_summary = cur;
}

int PhotonClusterBuilder::End(PHCompositeNode* /*topNode*/)
{
  if (m_iso_audit_mode)
  {
    print_audit_summary(true);
  }
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
    static unsigned long long s_evt = 0;
    ++s_evt;
    ++m_evt_seen;

    auto finish_event = [&](int rc)
    {
      if (m_iso_audit_mode)
      {
        print_audit_summary(false);
      }
      return rc;
    };

    if (!m_rawclusters)
    {
      m_rawclusters = findNode::getClass<RawClusterContainer>(topNode, m_input_cluster_node);
      if (!m_rawclusters)
      {
        ++m_evt_missing_rawclusters;
        std::cerr << Name() << ": missing RawClusterContainer '" << m_input_cluster_node << "'" << std::endl;
        return finish_event(Fun4AllReturnCodes::ABORTEVENT);
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
    const float vzCutForInfo = m_vz_cut_cm;

    m_vertex = std::numeric_limits<float>::quiet_NaN();
    const char* vtx_source = "NONE";
    bool used_mbd_vertex = false;
    bool used_global_vertex = false;

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
      used_mbd_vertex = true;
    }

    // 2) Optional fallback: GlobalVertexMap (still reco, only if finite)
    if (!std::isfinite(m_vertex) && std::isfinite(gv_z))
    {
      m_vertex = gv_z;
      vtx_source = "GlobalVertex";
      used_global_vertex = true;
    }

    // If still invalid, skip
    if (!std::isfinite(m_vertex))
    {
      ++m_evt_no_finite_vertex;
      if (Verbosity() >= 1)
      {
        std::cout << Name()
                  << ": [evt=" << s_evt << "] SKIP (no finite reco vertex)"
                  << " | MbdVertexMap n=" << mbd_n << " z=" << (std::isfinite(mbd_z) ? std::to_string(mbd_z) : std::string("NaN/NA"))
                  << " | GlobalVertexMap n=" << gv_n << " z=" << (std::isfinite(gv_z) ? std::to_string(gv_z) : std::string("NaN/NA"))
                  << std::endl;
      }
      return finish_event(Fun4AllReturnCodes::EVENT_OK);
    }

    if (used_mbd_vertex) ++m_evt_vertex_from_mbd;
    else if (used_global_vertex) ++m_evt_vertex_from_global;

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
      ++m_evt_skip_vz;
      if (Verbosity() >= 1)
      {
        std::cout << Name()
                  << ": [evt=" << s_evt << "] SKIP (|vz| cut)"
                  << " | vz=" << m_vertex << " | cut=" << vzCutForInfo
                  << " | source=" << vtx_source
                  << std::endl;
      }
      return finish_event(Fun4AllReturnCodes::EVENT_OK);
    }

  // iterate over clusters via map to have access to keys if needed
  const auto& rcmap = m_rawclusters->getClustersMap();

    size_t nClusters = rcmap.size();
    size_t nPassET = 0;
    size_t nBuilt = 0;
    size_t nMeanTimeDefault = 0;
    size_t nNonFiniteEta = 0;
    size_t nNonFinitePhi = 0;
    size_t nNonFiniteET = 0;

    float bestE = -999.0f;
    float bestET = -999.0f;
    float bestEta = -999.0f;
    float bestPhi = -999.0f;
    int bestLeadIeta = -999;
    int bestLeadIphi = -999;
    std::size_t bestTowermapSize = 0;
    float best_et1 = -999.0f;
    float best_et2 = -999.0f;
    float best_et3 = -999.0f;
    float best_et4 = -999.0f;
    float best_avg_eta = -999.0f;
    float best_avg_phi = -999.0f;

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

      if (!std::isfinite(eta)) ++nNonFiniteEta;
      if (!std::isfinite(phi)) ++nNonFinitePhi;
      if (!std::isfinite(ET)) ++nNonFiniteET;

      if (std::isfinite(ET) && ET > bestET)
      {
        bestE = E;
        bestET = ET;
        bestEta = eta;
        bestPhi = phi;

        const auto leadtowerindex = rc->get_lead_tower();
        bestLeadIeta = leadtowerindex.first;
        bestLeadIphi = leadtowerindex.second;
        bestTowermapSize = rc->get_towermap().size();

        const std::vector<float> showershape_dbg = rc->get_shower_shapes(m_shape_min_tower_E);
        if (!showershape_dbg.empty())
        {
          if (showershape_dbg.size() > 0) best_et1 = showershape_dbg[0];
          if (showershape_dbg.size() > 1) best_et2 = showershape_dbg[1];
          if (showershape_dbg.size() > 2) best_et3 = showershape_dbg[2];
          if (showershape_dbg.size() > 3) best_et4 = showershape_dbg[3];
          if (showershape_dbg.size() > 4) best_avg_eta = showershape_dbg[4] + 0.5F;
          if (showershape_dbg.size() > 5) best_avg_phi = showershape_dbg[5] + 0.5F;
        }
        else
        {
          best_et1 = -999.0f;
          best_et2 = -999.0f;
          best_et3 = -999.0f;
          best_et4 = -999.0f;
          best_avg_eta = -999.0f;
          best_avg_phi = -999.0f;
        }
      }

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

      if (nBuilt == 0)
      {
        std::cout << Name()
                  << ": [evt=" << s_evt << "] ZERO-PHOTON DEBUG"
                  << " | nonFinite(eta,phi,ET)=(" << nNonFiniteEta
                  << "," << nNonFinitePhi
                  << "," << nNonFiniteET << ")";

        if (nClusters == 0)
        {
          std::cout << " | reason=no input clusters in " << m_input_cluster_node;
        }
        else if (nPassET == 0)
        {
          std::cout << " | likely reason=all clusters failed ET threshold"
                    << " | bestCluster: E=" << bestE
                    << " ET=" << bestET
                    << " eta=" << bestEta
                    << " phi=" << bestPhi
                    << " lead=(" << bestLeadIeta << "," << bestLeadIphi << ")"
                    << " towermap=" << bestTowermapSize
                    << " | shower(et1,et2,et3,et4)=("
                    << best_et1 << "," << best_et2 << "," << best_et3 << "," << best_et4 << ")"
                    << " | avgTower=(" << best_avg_eta << "," << best_avg_phi << ")";
        }
        else
        {
          std::cout << " | reason=clusters passed ET but no photons survived downstream shaping/storage"
                    << " | bestCluster: E=" << bestE
                    << " ET=" << bestET
                    << " eta=" << bestEta
                    << " phi=" << bestPhi
                    << " lead=(" << bestLeadIeta << "," << bestLeadIphi << ")"
                    << " towermap=" << bestTowermapSize
                    << " | shower(et1,et2,et3,et4)=("
                    << best_et1 << "," << best_et2 << "," << best_et3 << "," << best_et4 << ")"
                    << " | avgTower=(" << best_avg_eta << "," << best_avg_phi << ")";
        }

        std::cout << std::endl;
      }
    }

    m_clusters_seen_total += nClusters;
    m_clusters_pass_et_total += nPassET;
    m_photons_built_total += nBuilt;
    m_mean_time_default_total += nMeanTimeDefault;
    m_nonfinite_eta_total += nNonFiniteEta;
    m_nonfinite_phi_total += nNonFinitePhi;
    m_nonfinite_et_total += nNonFiniteET;

    if (nClusters == 0)
    {
      ++m_evt_zero_input_clusters;
    }
    else if (nPassET == 0)
    {
      ++m_evt_zero_pass_et_clusters;
    }
    else if (nBuilt == 0)
    {
      ++m_evt_zero_built_after_pass_et;
    }
    else
    {
      ++m_evt_built_photons;
    }

  return finish_event(Fun4AllReturnCodes::EVENT_OK);
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

  const bool use_variant_a_iso = (m_emc_tower_node.find("_PHOSUB") != std::string::npos);

  auto calculate_layer_iso_signed = [&](TowerInfoContainer* towerContainer,
                                          RawTowerGeomContainer* geomContainer,
                                          RawTowerDefs::CalorimeterId calo_id,
                                          float radius,
                                          float core_radius)
    {
      if (!towerContainer || !geomContainer)
      {
        return std::numeric_limits<float>::quiet_NaN();
      }

      float layer_et = 0.0f;
      const unsigned int ntowers = towerContainer->size();
      for (unsigned int channel = 0; channel < ntowers; ++channel)
      {
        TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
        if (!tower || !tower->get_isGood())
        {
          continue;
        }

        const unsigned int towerkey = towerContainer->encode_key(channel);
        const int ieta = towerContainer->getTowerEtaBin(towerkey);
        const int iphi = towerContainer->getTowerPhiBin(towerkey);

        const RawTowerDefs::keytype geom_key = RawTowerDefs::encode_towerid(calo_id, ieta, iphi);
        RawTowerGeom* tower_geom = geomContainer->get_tower_geometry(geom_key);
        if (!tower_geom)
        {
          continue;
        }

        const double tower_eta = getTowerEta(tower_geom, 0, 0, m_vertex);
        const double tower_phi = tower_geom->get_phi();
        if (!std::isfinite(tower_eta) || !std::isfinite(tower_phi))
        {
          continue;
        }

        const double dr = deltaR(cluster_eta, cluster_phi, tower_eta, tower_phi);
        if (dr >= radius)
        {
          continue;
        }
        if (core_radius > 0.0f && dr < core_radius)
        {
          continue;
        }

        const float energy = tower->get_energy();
        if (!std::isfinite(energy))
        {
          continue;
        }

        layer_et += energy / std::cosh(tower_eta);
      }

      return layer_et;
    };

    auto compute_layer_iso = [&](RawTowerDefs::CalorimeterId calo_id, float radius)
    {
      TowerInfoContainer* container = nullptr;
      RawTowerGeomContainer* geom = nullptr;
      RawTowerDefs::CalorimeterId geom_id = calo_id;

      if (calo_id == RawTowerDefs::CalorimeterId::CEMC)
      {
        if (!use_variant_a_iso && m_is_auau && m_emc_tower_container_iso && m_geomEM_iso)
        {
          container = m_emc_tower_container_iso;
          geom = m_geomEM_iso;
          geom_id = RawTowerDefs::CalorimeterId::HCALIN;
        }
        else
        {
          container = m_emc_tower_container;
          geom = m_geomEM;
        }
      }
      else if (calo_id == RawTowerDefs::CalorimeterId::HCALIN)
      {
        if (!use_variant_a_iso && m_is_auau && m_ihcal_tower_container_iso)
        {
          container = m_ihcal_tower_container_iso;
        }
        else
        {
          container = m_ihcal_tower_container;
        }
        geom = m_geomIH;
      }
      else
      {
        if (!use_variant_a_iso && m_is_auau && m_ohcal_tower_container_iso)
        {
          container = m_ohcal_tower_container_iso;
        }
        else
        {
          container = m_ohcal_tower_container;
        }
        geom = m_geomOH;
      }

      if (use_variant_a_iso)
      {
        const float core_radius =
            (calo_id == RawTowerDefs::CalorimeterId::CEMC && radius > 0.05f) ? 0.05f : 0.0f;
        return calculate_layer_iso_signed(container, geom, geom_id, radius, core_radius);
      }

      return calculate_layer_et(cluster_eta, cluster_phi, radius, container, geom, geom_id, m_vertex);
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
    // Build iso values
    // -----------------------------
    const float iso_04_emcal  = use_variant_a_iso ? emcal_et_04  : (emcal_et_04  - ET);
    const float iso_03_emcal  = use_variant_a_iso ? emcal_et_03  : (emcal_et_03  - ET);
    const float iso_02_emcal  = use_variant_a_iso ? emcal_et_02  : (emcal_et_02  - ET);
    const float iso_01_emcal  = use_variant_a_iso ? emcal_et_01  : (emcal_et_01  - ET);
    const float iso_005_emcal = use_variant_a_iso ? emcal_et_005 : (emcal_et_005 - ET);

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

              if (!std::isfinite(e))
              {
                  ++nBelowThr;
                  droppedInCone.push_back({ieta, iphi, dRval, e, 0.0, isGood, "nonfiniteE"});
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
                      << " | inR_nonFiniteE=" << nBelowThr
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

  unsigned int nNull = 0;
  unsigned int nIsGood = 0;
  unsigned int nIsBad = 0;
  unsigned int nGeomMissing = 0;
  unsigned int nInConeGood = 0;
  unsigned int nBelowThreshold = 0;
  unsigned int nAccepted = 0;

  const char* calo_name = "UNKNOWN";
  if (calo_id == RawTowerDefs::CalorimeterId::CEMC)
  {
    calo_name = "CEMC";
  }
  else if (calo_id == RawTowerDefs::CalorimeterId::HCALIN)
  {
    calo_name = "HCALIN";
  }
  else if (calo_id == RawTowerDefs::CalorimeterId::HCALOUT)
  {
    calo_name = "HCALOUT";
  }

  for (unsigned int channel = 0; channel < ntowers; ++channel)
  {
    TowerInfo* tower = towerContainer->get_tower_at_channel(channel);
    if (!tower)
    {
      ++nNull;
      continue;
    }

    if (tower->get_isGood())
    {
      ++nIsGood;
    }
    else
    {
      ++nIsBad;
      continue;
    }

    unsigned int towerkey = towerContainer->encode_key(channel);
    int ieta = towerContainer->getTowerEtaBin(towerkey);
    int iphi = towerContainer->getTowerPhiBin(towerkey);

    RawTowerDefs::keytype geom_key = RawTowerDefs::encode_towerid(calo_id, ieta, iphi);
    RawTowerGeom* tower_geom = geomContainer->get_tower_geometry(geom_key);
    if (!tower_geom)
    {
      ++nGeomMissing;
      continue;
    }

      double tower_eta = getTowerEta(tower_geom, 0, 0, vertex_z);
      double tower_phi = tower_geom->get_phi();
      if (!std::isfinite(tower_eta) || !std::isfinite(tower_phi))
      {
        continue;
      }

      if (deltaR(seed_eta, seed_phi, tower_eta, tower_phi) >= radius)
      {
        continue;
      }

      ++nInConeGood;

      float energy = tower->get_energy();
      if (!std::isfinite(energy))
      {
        ++nBelowThreshold;
        continue;
      }

      float et = energy / std::cosh(tower_eta);
      layer_et += et;
      ++nAccepted;
  }

  if (Verbosity() >= 4)
  {
    std::cout << Name()
              << ": calculate_layer_et isGood summary"
              << " | calo=" << calo_name
              << " | R=" << radius
              << " | ntowers=" << ntowers
              << " | null=" << nNull
              << " | isGood=1=" << nIsGood
              << " | isGood=0=" << nIsBad
              << " | geomMissing=" << nGeomMissing
              << " | inCone_good=" << nInConeGood
              << " | inCone_nonFiniteE=" << nBelowThreshold
              << " | accepted=" << nAccepted
              << " | sumEt=" << layer_et
              << std::endl;
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
