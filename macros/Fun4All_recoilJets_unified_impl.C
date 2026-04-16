//======================================================================
//  Fun4All_recoilJets_unified_impl.C
//  --------------------------------------------------------------------
#pragma once
#if !defined(RJ_UNIFIED_ANALYSIS_PP) && !defined(RJ_UNIFIED_ANALYSIS_AUAU)
  #error "Define RJ_UNIFIED_ANALYSIS_PP or RJ_UNIFIED_ANALYSIS_AUAU before including Fun4All_recoilJets_unified_impl.C"
#endif
#if defined(__CINT__) || defined(__CLING__)
  R__ADD_INCLUDE_PATH(/sphenix/u/patsfan753/thesisAnalysis/install/include)
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
  R__ADD_INCLUDE_PATH(/sphenix/u/patsfan753/thesisAnalysis_auau/install/include)
#endif
#endif
#if defined(__CLING__)
  #pragma cling add_include_path("/sphenix/u/patsfan753/thesisAnalysis/install/include")
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
  #pragma cling add_include_path("/sphenix/u/patsfan753/thesisAnalysis_auau/install/include")
#endif
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

//–––– Standard Fun4All ––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllNoSyncDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfo.h>
#include <jetbase/JetContainer.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <mbd/MbdPmtContainer.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerBuilder.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/CDBInterface.h>
#include <clusteriso/ClusterIso.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <calobase/RawTowerGeom.h>
#include <caloreco/RawTowerCalibration.h>
#include "/sphenix/u/patsfan753/thesisAnalysis/install/include/caloreco/PhotonClusterBuilder.h"
#include <jetbase/Jet.h>
#include <g4jets/TruthJetInput.h>

#include <jetbase/FastJetOptions.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <globalvertex/GlobalVertexReco.h>
#include <caloreco/CaloTowerCalib.h>

#include <centrality/CentralityReco.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <jetbase/FastJetOptions.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbase/JetCalib.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/TowerBackground.h>
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
#include <eventplaneinfo/EventPlaneReco.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src_AuAu/RecoilJets_AuAu.h"
#else
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h"
#endif

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>     // getenv
#include <algorithm>   // std::transform
#include <cctype>      // std::tolower
#include <iomanip>     // std::setw, std::setprecision
#include <dlfcn.h>     // dlopen RTLD_NOLOAD
#include <typeinfo>    // typeid RTTI probe
#include <TSystem.h>   // gSystem, GetBuildArch/Compiler info
#include <csignal>     // signal handlers (debug backtrace)
#include <execinfo.h>  // backtrace
#include <unistd.h>    // STDERR_FILENO
#include <cstdio>      // snprintf
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/macros/Calo_Calib.C"

// Load local CaloReco/CaloIO first so PhotonClusterBuilder and CaloReco stay on your thesisAnalysis install.
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis_auau/install/lib/libRecoilJetsAuAu.so)
#else
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)
#endif
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)


// Then load the rest of the environment stack
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libglobalvertex.so)
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
R__LOAD_LIBRARY(libeventplaneinfo.so)
#endif
R__LOAD_LIBRARY(libcentrality.so)      // always
R__LOAD_LIBRARY(libcentrality_io.so)   // if you instantiate CentralityReco
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libzdcinfo.so)
R__LOAD_LIBRARY(libmbd.so)

//======================================================================
//  Convenience helpers
//======================================================================
namespace detail
{
  /// Throw a nicely formatted exception on unrecoverable error
  [[noreturn]] void bail(const std::string& msg)
  {
    std::ostringstream oss;
    oss << "\n[FATAL] Fun4All_recoilJets :: " << msg << '\n';
    throw std::runtime_error(oss.str());
  }

  /// Trim whitespace from both ends (for robust list-file parsing)
  inline std::string trim(std::string s)
  {
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
  }
}

namespace detail
{
  inline FastJetAlgoSub* fjAlgo(const float R)
  {
    FastJetOptions o{};               // IMPORTANT: value-initialize ALL fields to safe defaults
    o.algo            = Jet::ANTIKT;  // algorithm
    o.jet_R           = R;            // jet radius
    o.use_jet_min_pt  = true;         // enable a ptmin
    o.jet_min_pt      = 0.0f;         // ptmin value
    o.verbosity       = 0;            // quiet
    return new FastJetAlgoSub(o);
  }
}

namespace yamlcfg
{
  // Forward declarations (LoadJetRKeys appears before helper definitions below)
  inline std::string ResolveYAMLPath();
  inline bool ReadWholeFile(const std::string& path, std::string& out);
  inline bool StartsWithKey(const std::string& line, const std::string& key);
  inline std::string AfterColon(const std::string& line);
  inline void ParseInlineListDoubles(std::string s, std::vector<double>& out);

  inline std::vector<std::string> LoadJetRKeys(int vlevel)
  {
    std::vector<std::string> outKeys;

    const std::string yamlPath = ResolveYAMLPath();
    std::string yamlText;
    const bool ok = ReadWholeFile(yamlPath, yamlText);

    // Default behavior (baseline): r02 + r04
    std::vector<double> radii = {0.2, 0.4};

    if (ok)
    {
      std::istringstream iss(yamlText);
      for (std::string line; std::getline(iss, line); )
      {
        line = detail::trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;

        if (StartsWithKey(line, "jet_radii"))
        {
          const std::string rhs = AfterColon(line);
          ParseInlineListDoubles(rhs, radii);
          break;
        }
      }
    }

    auto push_unique = [&](const std::string& k)
    {
      if (std::find(outKeys.begin(), outKeys.end(), k) == outKeys.end())
        outKeys.push_back(k);
    };

    for (double R : radii)
    {
      if (!std::isfinite(R) || R <= 0.0) continue;

      // Map R to "r%02d" where D ~ round(10*R)
      const int D = static_cast<int>(R * 10.0 + 0.5);
      if (D <= 0) continue;

      std::ostringstream rk;
      rk << "r" << std::setw(2) << std::setfill('0') << D;
      push_unique(rk.str());
    }

    if (outKeys.empty())
    {
      outKeys = {"r02","r04"};
      if (vlevel > 0)
        std::cout << "[CFG] jet_radii: no valid entries found -> defaulting to [r02,r04]\n";
    }

    if (vlevel > 0)
    {
      std::cout << "[CFG] jet_radii from YAML (" << yamlPath << "): [";
      for (std::size_t i = 0; i < outKeys.size(); ++i)
        std::cout << outKeys[i] << (i + 1 < outKeys.size() ? ", " : "");
      std::cout << "]\n";
    }

    return outKeys;
  }

  inline bool WantRKey(const std::vector<std::string>& keys, const std::string& key)
  {
    if (keys.empty()) return true;
    return (std::find(keys.begin(), keys.end(), key) != keys.end());
  }
}

namespace yamlcfg
{
    struct Config
    {
        // file provenance
        std::string yamlPath = "";
        std::string yamlText = "";
        
        // baseline defaults (match current behavior)
        double photon_eta_abs_max = 0.7;
        double jet_pt_min = 10.0;
        double back_to_back_dphi_min_pi_fraction = 0.875;
        
        bool   use_vz_cut = true;
        double vz_cut_cm  = 30.0;
        
        std::vector<int> centrality_edges = {0, 10, 20, 40, 60, 80, 100};
        
        bool vertex_reweight_on = false;
        std::string vertex_reweight_file = "";
        std::string vertex_reweight_hist = "data_over_MC_ratios/h_zvtx_ratio_data_over_photonJet";
        
        bool centrality_reweight_on = false;
        std::string centrality_reweight_file = "/sphenix/u/patsfan753/scratch/thesisAnalysis/reweightingDer/output/centrality_reweighting.root";
        std::string centrality_reweight_hist = "nom_cent_rw_hist";
        
        double isoA = 1.08128;
        double isoB = 0.0299107;
        double isoGap = 1.0;
        double isoFixed = 2.0;
        double truthIsoGeV = 4.0;
        double isoConeR = 0.30;
        double isoTowMin = 0.0;
        bool   isSlidingIso = true;
        bool   isSlidingAndFixed = false;
        
        // Per-centrality AuAu sliding WPs: each entry = {aGeV, bPerGeV, sideGapGeV}
        struct CentIsoWP { double aGeV; double bPerGeV; double sideGapGeV; };
        std::vector<CentIsoWP> auauCentIsoWP;
        
        // Photon ID cuts (PPG12 Table 4) baseline
        double pre_e11e33_max = 0.98;
        double pre_et1_min    = 0.60;
        double pre_et1_max    = 1.00;
        double pre_e32e35_min = 0.80;
        double pre_e32e35_max = 1.00;
        double pre_weta_max   = 0.60;
        
        double tight_w_lo            = 0.0;
        double tight_w_hi_intercept  = 0.15;
        double tight_w_hi_slope      = 0.006;
        
        double tight_e11e33_min = 0.40;
        double tight_e11e33_max = 0.98;
        
        double tight_et1_min    = 0.90;
        double tight_et1_max    = 1.00;
        
        double tight_e32e35_min = 0.92;
        double tight_e32e35_max = 1.00;
        
        double pho_dr_max = 0.05;
        double jet_dr_max = 0.3;
        
        std::vector<double> jes3_photon_pt_bins = {15,17,19,21,23,26,35};
        std::vector<double> unfold_reco_photon_pt_bins  = {10,15,17,19,21,23,26,35,40};
        std::vector<double> unfold_truth_photon_pt_bins = {5,10,15,17,19,21,23,26,35,40};
        
        double unfold_jet_pt_start = 0.0;
        double unfold_jet_pt_stop  = 60.0;
        double unfold_jet_pt_step  = 0.5;
        
        std::vector<double> unfold_xj_bins = {0.0,0.20,0.24,0.29,0.35,0.41,0.50,0.60,0.72,0.86,1.03,1.24,1.49,1.78,2.14,3.0};
        
        // EventDisplay diagnostics payload (EventDisplayTree)
        bool event_display_tree = true;
        int  event_display_tree_max_per_bin = 0;
        std::string clusterUEpipeline = "noSub";
        bool doPi0Analysis = false;
    };
    
    inline std::string DefaultYAMLPath()
    {
        std::string here = __FILE__;
        const std::size_t p = here.find_last_of('/');
        if (p == std::string::npos) return std::string("analysis_config.yaml");
        return here.substr(0, p) + "/analysis_config.yaml";
    }
    
    inline std::string ResolveYAMLPath()
    {
        if (const char* env = std::getenv("RJ_CONFIG_YAML"))
        {
            std::string s = detail::trim(std::string(env));
            if (!s.empty()) return s;
        }
        return DefaultYAMLPath();
    }
    
    inline bool ReadWholeFile(const std::string& path, std::string& out)
    {
        std::ifstream in(path);
        if (!in.is_open()) return false;
        std::ostringstream ss;
        ss << in.rdbuf();
        out = ss.str();
        return true;
    }
    
    inline bool StartsWithKey(const std::string& line, const std::string& key)
    {
        const std::string k = key + ":";
        if (line.size() < k.size()) return false;
        return (line.compare(0, k.size(), k) == 0);
    }
    
    inline std::string AfterColon(const std::string& line)
    {
        const std::size_t c = line.find(':');
        if (c == std::string::npos) return std::string{};
        return detail::trim(line.substr(c + 1));
    }
    
    inline bool ParseBool(std::string s, bool& out)
    {
        s = detail::trim(s);
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        if (s == "true"  || s == "1" || s == "yes" || s == "on")  { out = true;  return true; }
        if (s == "false" || s == "0" || s == "no"  || s == "off") { out = false; return true; }
        return false;
    }
    
    inline bool ParseDouble(const std::string& s, double& out)
    {
        try { out = std::stod(detail::trim(s)); return true; } catch (...) { return false; }
    }
    
    inline bool ParseInt(const std::string& s, int& out)
    {
        try { out = std::stoi(detail::trim(s)); return true; } catch (...) { return false; }
    }
    
    inline void ParseInlineListDoubles(std::string s, std::vector<double>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('[');
        const std::size_t r = s.find(']');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            double v = 0.0;
            if (ParseDouble(tok, v)) out.push_back(v);
        }
    }
    
    inline void ParseInlineListInts(std::string s, std::vector<int>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('[');
        const std::size_t r = s.find(']');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string tok;
        while (std::getline(ss, tok, ','))
        {
            int v = 0;
            if (ParseInt(tok, v)) out.push_back(v);
        }
    }
    
    inline void ParseInlineMapDoubles(std::string s, std::map<std::string, double>& out)
    {
        out.clear();
        s = detail::trim(s);
        const std::size_t l = s.find('{');
        const std::size_t r = s.find('}');
        if (l == std::string::npos || r == std::string::npos || r <= l) return;
        std::string inner = s.substr(l + 1, r - l - 1);
        std::stringstream ss(inner);
        std::string pair;
        while (std::getline(ss, pair, ','))
        {
            const std::size_t c = pair.find(':');
            if (c == std::string::npos) continue;
            std::string k = detail::trim(pair.substr(0, c));
            std::string v = detail::trim(pair.substr(c + 1));
            double dv = 0.0;
            if (ParseDouble(v, dv)) out[k] = dv;
        }
    }
    
    inline Config LoadConfig()
    {
        Config cfg;
        cfg.yamlPath = ResolveYAMLPath();
        
        const bool okRead = ReadWholeFile(cfg.yamlPath, cfg.yamlText);
        
        // Local verbosity for YAML parsing diagnostics (independent of the later banner)
        int yamlV = 0;
        if (const char* venv = std::getenv("RJ_VERBOSITY"))
        {
            try { yamlV = std::stoi(venv); } catch (...) { yamlV = 0; }
        }
        
        bool strict = false;
        if (const char* env = std::getenv("RJ_YAML_STRICT")) strict = (std::atoi(env) != 0);
        
        auto warn_parse = [&](const std::string& key, const std::string& rhs, const std::string& why)
        {
            if (yamlV > 0)
            {
                std::cout << "[CFG][WARN] YAML parse failed for key '" << key << "'"
                << " (rhs='" << rhs << "')"
                << " -> " << why
                << " (keeping default)\n";
            }
            if (strict)
            {
                std::ostringstream oss;
                oss << "YAML parse failed for key '" << key << "' (rhs='" << rhs << "'): " << why;
                throw std::runtime_error(oss.str());
            }
        };
        
        auto info_parse = [&](const std::string& msg)
        {
            if (yamlV > 0) std::cout << msg << "\n";
        };
        
        if (!okRead || cfg.yamlText.empty())
        {
            if (yamlV > 0)
            {
                std::cout << "[CFG][WARN] Could not read YAML (or file empty): " << cfg.yamlPath
                << " (all values will remain defaults)\n";
            }
            if (strict)
            {
                throw std::runtime_error(std::string("Could not read YAML: ") + cfg.yamlPath);
            }
            return cfg;
        }
        
        std::istringstream iss(cfg.yamlText);
        for (std::string line; std::getline(iss, line); )
        {
            line = detail::trim(line);
            if (line.empty()) continue;
            if (!line.empty() && line[0] == '#') continue;
            
            if (StartsWithKey(line, "photon_eta_abs_max"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.photon_eta_abs_max))
                    warn_parse("photon_eta_abs_max", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "jet_pt_min"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.jet_pt_min))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.jet_pt_min = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] jet_pt_min is a list (n=" << v.size() << "); using first value = " << cfg.jet_pt_min;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("jet_pt_min", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "back_to_back_dphi_min_pi_fraction"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.back_to_back_dphi_min_pi_fraction))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.back_to_back_dphi_min_pi_fraction = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] back_to_back_dphi_min_pi_fraction is a list (n=" << v.size()
                        << "); using first value = " << cfg.back_to_back_dphi_min_pi_fraction;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("back_to_back_dphi_min_pi_fraction", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "use_vz_cut"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.use_vz_cut))
                    warn_parse("use_vz_cut", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "vz_cut_cm"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.vz_cut_cm))
                    warn_parse("vz_cut_cm", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "coneR"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.isoConeR))
                    warn_parse("coneR", rhs, "expected a scalar double");
            }
            else if (StartsWithKey(line, "centrality_edges"))
            {
              const std::string rhs = AfterColon(line);
              std::vector<int> v;
              ParseInlineListInts(rhs, v);
              if (v.size() >= 2)
              {
                cfg.centrality_edges = v;
              }
              else
              {
                warn_parse("centrality_edges", rhs, "expected an inline list of integers with size >= 2 (e.g. [0, 10, 20, 40, 60, 80, 100])");
              }
            }
            else if (StartsWithKey(line, "vertex_reweight_on"))
            {
              const std::string rhs = AfterColon(line);
              if (!ParseBool(rhs, cfg.vertex_reweight_on))
                warn_parse("vertex_reweight_on", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "vertex_reweight_file"))
            {
              cfg.vertex_reweight_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "vertex_reweight_hist"))
            {
              cfg.vertex_reweight_hist = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "centrality_reweight_on"))
            {
              const std::string rhs = AfterColon(line);
              if (!ParseBool(rhs, cfg.centrality_reweight_on))
                warn_parse("centrality_reweight_on", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "centrality_reweight_file"))
            {
              cfg.centrality_reweight_file = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "centrality_reweight_hist"))
            {
              cfg.centrality_reweight_hist = detail::trim(AfterColon(line));
            }
            else if (StartsWithKey(line, "clusterUEpipeline"))
            {
                std::string rhs = AfterColon(line);
                // trim leading/trailing whitespace
                while (!rhs.empty() && (rhs.front() == ' ' || rhs.front() == '\t')) rhs.erase(rhs.begin());
                while (!rhs.empty() && (rhs.back() == ' ' || rhs.back() == '\t')) rhs.pop_back();
                if (rhs == "true" || rhs == "1")
                    cfg.clusterUEpipeline = "variantA";
                else if (rhs == "false" || rhs == "0")
                    cfg.clusterUEpipeline = "noSub";
                else if (rhs == "noSub" || rhs == "baseVariant" || rhs == "variantA" || rhs == "variantB")
                    cfg.clusterUEpipeline = rhs;
                else
                    warn_parse("clusterUEpipeline", rhs, "expected noSub|baseVariant|variantA|variantB (or true/false for compat)");
            }
            else if (StartsWithKey(line, "doPi0Analysis"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.doPi0Analysis))
                    warn_parse("doPi0Analysis", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "event_display_tree"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.event_display_tree))
                    warn_parse("event_display_tree", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "event_display_tree_max_per_bin"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseInt(rhs, cfg.event_display_tree_max_per_bin))
                {
                    warn_parse("event_display_tree_max_per_bin", rhs, "expected an integer");
                }
                if (cfg.event_display_tree_max_per_bin < 0) cfg.event_display_tree_max_per_bin = 0;
            }
            else if (StartsWithKey(line, "matching"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("pho_dr_max")) cfg.pho_dr_max = m["pho_dr_max"];
                if (m.count("jet_dr_max")) cfg.jet_dr_max = m["jet_dr_max"];
            }
            else if (StartsWithKey(line, "isolation_wp"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("aGeV"))       cfg.isoA        = m["aGeV"];
                if (m.count("bPerGeV"))    cfg.isoB        = m["bPerGeV"];
                if (m.count("sideGapGeV")) cfg.isoGap      = m["sideGapGeV"];
                if (m.count("fixedGeV"))   cfg.isoFixed    = m["fixedGeV"];
                if (m.count("truthIsoGeV"))cfg.truthIsoGeV = m["truthIsoGeV"];
                if (m.count("coneR"))      cfg.isoConeR    = m["coneR"];
                if (m.count("towerMin"))   cfg.isoTowMin   = m["towerMin"];
            }
            else if (StartsWithKey(line, "isSlidingIso"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.isSlidingIso))
                    warn_parse("isSlidingIso", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "isSlidingAndFixed"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseBool(rhs, cfg.isSlidingAndFixed))
                    warn_parse("isSlidingAndFixed", rhs, "expected true/false");
            }
            else if (StartsWithKey(line, "fixedGeV"))
            {
                const std::string rhs = AfterColon(line);
                if (!ParseDouble(rhs, cfg.isoFixed))
                {
                    std::vector<double> v;
                    ParseInlineListDoubles(rhs, v);
                    if (!v.empty())
                    {
                        cfg.isoFixed = v.front();
                        std::ostringstream oss;
                        oss << "[CFG] fixedGeV is a list (n=" << v.size() << "); using first value = " << cfg.isoFixed;
                        info_parse(oss.str());
                    }
                    else
                    {
                        warn_parse("fixedGeV", rhs, "expected a scalar double or an inline list [..]");
                    }
                }
            }
            else if (StartsWithKey(line, "auau_cent_iso_wp"))
            {
                // Multi-line list: read subsequent "  - {aGeV: ..., bPerGeV: ..., sideGapGeV: ...}" lines
                cfg.auauCentIsoWP.clear();
                while (std::getline(iss, line))
                {
                    line = detail::trim(line);
                    if (line.empty()) continue;
                    if (line[0] == '#') continue;
                    if (line[0] != '-') break;  // end of list
                    std::map<std::string, double> m;
                    ParseInlineMapDoubles(line.substr(1), m);
                    Config::CentIsoWP wp{};
                    wp.aGeV       = m.count("aGeV")       ? m["aGeV"]       : cfg.isoA;
                    wp.bPerGeV    = m.count("bPerGeV")     ? m["bPerGeV"]    : cfg.isoB;
                    wp.sideGapGeV = m.count("sideGapGeV")  ? m["sideGapGeV"] : cfg.isoGap;
                    cfg.auauCentIsoWP.push_back(wp);
                }
                if (yamlV > 0)
                {
                    std::cout << "[CFG] auau_cent_iso_wp: parsed " << cfg.auauCentIsoWP.size() << " entries\n";
                }
            }
            else if (StartsWithKey(line, "photon_id_pre"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("e11e33_max")) cfg.pre_e11e33_max = m["e11e33_max"];
                if (m.count("et1_min"))    cfg.pre_et1_min    = m["et1_min"];
                if (m.count("et1_max"))    cfg.pre_et1_max    = m["et1_max"];
                if (m.count("e32e35_min")) cfg.pre_e32e35_min = m["e32e35_min"];
                if (m.count("e32e35_max")) cfg.pre_e32e35_max = m["e32e35_max"];
                if (m.count("weta_max"))   cfg.pre_weta_max   = m["weta_max"];
            }
            else if (StartsWithKey(line, "photon_id_tight"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("w_lo"))            cfg.tight_w_lo           = m["w_lo"];
                if (m.count("w_hi_intercept"))  cfg.tight_w_hi_intercept = m["w_hi_intercept"];
                if (m.count("w_hi_slope"))      cfg.tight_w_hi_slope     = m["w_hi_slope"];
                
                if (m.count("e11e33_min")) cfg.tight_e11e33_min = m["e11e33_min"];
                if (m.count("e11e33_max")) cfg.tight_e11e33_max = m["e11e33_max"];
                
                if (m.count("et1_min"))    cfg.tight_et1_min    = m["et1_min"];
                if (m.count("et1_max"))    cfg.tight_et1_max    = m["et1_max"];
                
                if (m.count("e32e35_min")) cfg.tight_e32e35_min = m["e32e35_min"];
                if (m.count("e32e35_max")) cfg.tight_e32e35_max = m["e32e35_max"];
            }
            else if (StartsWithKey(line, "jes3_photon_pt_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.jes3_photon_pt_bins);
                if (cfg.jes3_photon_pt_bins.size() < 2)
                    warn_parse("jes3_photon_pt_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_reco_photon_pt_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.unfold_reco_photon_pt_bins);
                if (cfg.unfold_reco_photon_pt_bins.size() < 2)
                    warn_parse("unfold_reco_photon_pt_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_truth_photon_pt_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.unfold_truth_photon_pt_bins);
                if (cfg.unfold_truth_photon_pt_bins.size() < 2)
                    warn_parse("unfold_truth_photon_pt_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_xj_bins"))
            {
                const std::string rhs = AfterColon(line);
                ParseInlineListDoubles(rhs, cfg.unfold_xj_bins);
                if (cfg.unfold_xj_bins.size() < 2)
                    warn_parse("unfold_xj_bins", rhs, "expected an inline list with >=2 edges");
            }
            else if (StartsWithKey(line, "unfold_jet_pt_binning"))
            {
                std::map<std::string, double> m;
                ParseInlineMapDoubles(AfterColon(line), m);
                if (m.count("start")) cfg.unfold_jet_pt_start = m["start"];
                if (m.count("stop"))  cfg.unfold_jet_pt_stop  = m["stop"];
                if (m.count("step"))  cfg.unfold_jet_pt_step  = m["step"];
            }
        }
        
        if (yamlV > 0)
        {
            std::cout << "[CFG] YAML parsing completed OK: " << cfg.yamlPath
            << (strict ? " (strict)" : "") << "\n";
        }
        
        return cfg;
    }
    
    inline void ExpandUniformEdges(std::vector<double>& edges, double start, double stop, double step)
    {
        edges.clear();
        if (!(std::isfinite(start) && std::isfinite(stop) && std::isfinite(step))) return;
        if (step <= 0.0) return;
        if (stop <= start) return;
        
        const int n = (int) std::llround((stop - start) / step);
        edges.reserve((std::size_t)n + 2);
        for (int i = 0; i <= n; ++i)
        {
            edges.push_back(start + step * (double)i);
        }
        if (edges.empty() || std::fabs(edges.back() - stop) > 1e-9) edges.push_back(stop);
    }
}

// ============================================================================
// Ensure the composite nodes JetCalib hard-requires exist.
// JetCalib::CreateNodeTree requires BOTH:
//   - PHCompositeNode "ANTIKT"
//   - PHCompositeNode "TOWER"
// Your pp DSTs often don't have "TOWER", so JetCalib aborts the run.
// ============================================================================
class EnsureJetCalibNodes final : public SubsysReco
{
 public:
  explicit EnsureJetCalibNodes(const std::string& name="EnsureJetCalibNodes")
    : SubsysReco(name) {}

  int InitRun(PHCompositeNode* topNode) override
  {
    PHNodeIterator iter(topNode);

    auto* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      std::cerr << "[EnsureJetCalibNodes] FATAL: DST node not found\n";
      return Fun4AllReturnCodes::ABORTRUN;
    }

    // JetCalib requires ANTIKT to exist
    auto* antiktNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "ANTIKT"));
    if (!antiktNode)
    {
      antiktNode = new PHCompositeNode("ANTIKT");
      dstNode->addNode(antiktNode);
      if (Verbosity() > 0)
        std::cout << "[EnsureJetCalibNodes] created DST/ANTIKT composite node\n";
    }

    // JetCalib requires TOWER to exist (this is what you're missing)
    auto* towerNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "TOWER"));
    if (!towerNode)
    {
      towerNode = new PHCompositeNode("TOWER");
      dstNode->addNode(towerNode);
      if (Verbosity() > 0)
        std::cout << "[EnsureJetCalibNodes] created DST/TOWER composite node\n";
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }
};

class ProcessEnvSetter final : public SubsysReco
{
 public:
  ProcessEnvSetter(const std::string& name,
                   const std::string& key,
                   const std::string& value)
    : SubsysReco(name)
    , m_key(key)
    , m_value(value)
  {}

  int InitRun(PHCompositeNode* /*topNode*/) override
  {
    setenv(m_key.c_str(), m_value.c_str(), 1);
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* /*topNode*/) override
  {
    setenv(m_key.c_str(), m_value.c_str(), 1);
    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_key;
  std::string m_value;
};


class TowerAudit final : public SubsysReco
{
 public:
  explicit TowerAudit(const std::string& name,
                      const std::string& nodeCEMC,
                      const std::string& nodeHCIN,
                      const std::string& nodeHCOUT,
                      int maxPrintEvents = 5)
    : SubsysReco(name)
    , m_nodeCEMC(nodeCEMC)
    , m_nodeHCIN(nodeHCIN)
    , m_nodeHCOUT(nodeHCOUT)
    , m_maxPrint(maxPrintEvents)
  {}

  int InitRun(PHCompositeNode* topNode) override
  {
    // Just verify nodes exist at InitRun (they may still be filled per-event)
    auto* cemc = findNode::getClass<TowerInfoContainer>(topNode, m_nodeCEMC);
    auto* hcin = findNode::getClass<TowerInfoContainer>(topNode, m_nodeHCIN);
    auto* hcot = findNode::getClass<TowerInfoContainer>(topNode, m_nodeHCOUT);

    if (Verbosity() > 0)
    {
      std::cout << "[" << Name() << "::InitRun]"
                << " nodes:"
                << " CEMC=" << m_nodeCEMC << (cemc ? "(OK)" : "(MISSING)")
                << " HCALIN=" << m_nodeHCIN << (hcin ? "(OK)" : "(MISSING)")
                << " HCALOUT=" << m_nodeHCOUT << (hcot ? "(OK)" : "(MISSING)")
                << std::endl;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  }

  int process_event(PHCompositeNode* topNode) override
  {
    ++m_evt;
    if (m_evt > m_maxPrint) return Fun4AllReturnCodes::EVENT_OK;

    auto auditOne = [&](const char* label, const std::string& node)
    {
      auto* cont = findNode::getClass<TowerInfoContainer>(topNode, node);
      if (!cont)
      {
        std::cout << "[" << Name() << "] evt=" << m_evt
                  << " " << label << " node=" << node << " MISSING\n";
        return;
      }

      const unsigned int nt = cont->size();
      unsigned int nNull = 0, nGood = 0, nBad = 0;
      unsigned int nEpos = 0, nEposGood = 0, nEposBad = 0;

      double sumE = 0.0, sumEgood = 0.0, sumEbad = 0.0;
      double sumT = 0.0, sumTgood = 0.0, sumTbad = 0.0;
      unsigned int nT = 0, nTgood = 0, nTbad = 0;

      // sample the whole container (24576 for CEMC is fine for a few events)
      for (unsigned int ch = 0; ch < nt; ++ch)
      {
        TowerInfo* t = cont->get_tower_at_channel(ch);
        if (!t) { ++nNull; continue; }

        const bool good = t->get_isGood();
        const float e   = t->get_energy();
        const float tim = t->get_time();

        if (good) ++nGood; else ++nBad;

        if (std::isfinite(e))
        {
          sumE += e;
          if (e > 0) ++nEpos;
          if (good) { sumEgood += e; if (e > 0) ++nEposGood; }
          else      { sumEbad  += e; if (e > 0) ++nEposBad;  }
        }

        if (std::isfinite(tim))
        {
          ++nT;
          sumT += tim;
          if (good) { ++nTgood; sumTgood += tim; }
          else      { ++nTbad;  sumTbad  += tim; }
        }
      }

      const double fracBad = (nt > 0) ? (double)nBad / (double)nt : 0.0;
      const double meanT   = (nT > 0) ? sumT / (double)nT : 0.0;
      const double meanTg  = (nTgood > 0) ? sumTgood / (double)nTgood : 0.0;
      const double meanTb  = (nTbad > 0) ? sumTbad / (double)nTbad : 0.0;

      std::cout << "[" << Name() << "] evt=" << m_evt << " " << label
                << " node=" << node
                << " ntowers=" << nt
                << " null=" << nNull
                << " good=" << nGood
                << " bad=" << nBad
                << " fracBad=" << std::fixed << std::setprecision(4) << fracBad
                << " | sumE=" << std::setprecision(3) << sumE
                << " (good=" << sumEgood << ", bad=" << sumEbad << ")"
                << " | E>0: all=" << nEpos << " good=" << nEposGood << " bad=" << nEposBad
                << " | meanTime=" << std::setprecision(3) << meanT
                << " (good=" << meanTg << ", bad=" << meanTb << ")"
                << "\n";
    };

    auditOne("CEMC",   m_nodeCEMC);
    auditOne("HCALIN", m_nodeHCIN);
    auditOne("HCALOUT",m_nodeHCOUT);

    return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_nodeCEMC;
  std::string m_nodeHCIN;
  std::string m_nodeHCOUT;
  int m_maxPrint = 5;
  int m_evt = 0;
};



class JetCalibOneEventProbe final : public SubsysReco
{
 public:
  JetCalibOneEventProbe(const std::string& name,
                        const std::string& rawNode,
                        const std::string& calibNode,
                        int maxJetsToPrint = 12)
    : SubsysReco(name)
    , m_rawNode(rawNode)
    , m_calibNode(calibNode)
    , m_maxJets(maxJetsToPrint)
  {}

  int process_event(PHCompositeNode* topNode) override
  {
    // Hard gate: ONLY show anything at ultra-verbose
    if (Verbosity() < 20) return Fun4AllReturnCodes::EVENT_OK;

    // Print only once for the whole job
    if (m_fired) return Fun4AllReturnCodes::EVENT_OK;

    auto* raw   = findNode::getClass<JetContainer>(topNode, m_rawNode);
    auto* calib = findNode::getClass<JetContainer>(topNode, m_calibNode);

      // Only trigger when both exist and actually have jets
      if (!raw || !calib) return Fun4AllReturnCodes::EVENT_OK;
      if (raw->size() == 0 || calib->size() == 0) return Fun4AllReturnCodes::EVENT_OK;

      // Grab z-vertex (MBD) JetCalib uses; REQUIRE |zvtx| < 30 cm
      float zvtx = -999.0f;
      auto* vtxmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
      if (vtxmap && !vtxmap->empty())
      {
        GlobalVertex* vtx = vtxmap->begin()->second;
        if (vtx)
        {
          auto mbdStart = vtx->find_vertexes(GlobalVertex::MBD);
          auto mbdEnd   = vtx->end_vertexes();
          for (auto it = mbdStart; it != mbdEnd; ++it)
          {
            const auto& [type, vec] = *it;
            if (type != GlobalVertex::MBD) continue;
            for (const auto* vv : vec)
            {
              if (!vv) continue;
              zvtx = vv->get_z();
            }
          }
        }
      }

      // Require a well-defined vertex AND a central-ish event
      if (!std::isfinite(zvtx)) return Fun4AllReturnCodes::EVENT_OK;
      if (std::fabs(zvtx) >= 30.0f) return Fun4AllReturnCodes::EVENT_OK;

      const int nRaw   = (int) raw->size();
      const int nCalib = (int) calib->size();
      const int nScan  = std::min(nRaw, nCalib);

      // Only print jets with pt_raw >= 5 GeV
      const float ptMinPrint = 5.0f;

      // First pass: count how many jets pass pt threshold (to decide whether to print at all)
      int nPass = 0;
      for (int i = 0; i < nScan; ++i)
      {
        const Jet* jr = raw->get_jet(i);
        if (!jr) continue;
        if (jr->get_pt() >= ptMinPrint) ++nPass;
      }

      // If nothing interesting, don't fire and keep searching future events
      if (nPass == 0) return Fun4AllReturnCodes::EVENT_OK;

      // ---------------- ANSI helpers ----------------
      const char* RST  = "\033[0m";
      const char* BOLD = "\033[1m";
      const char* DIM  = "\033[2m";

      const char* C_RAW   = "\033[38;5;45m";   // cyan-ish
      const char* C_CAL   = "\033[38;5;208m";  // orange-ish
      const char* C_INFO  = "\033[38;5;111m";  // light blue
      const char* C_WARN  = "\033[38;5;220m";  // yellow
      const char* C_GOOD  = "\033[38;5;82m";   // green
      const char* C_BAD   = "\033[38;5;196m";  // red
      const char* C_BAR   = "\033[38;5;244m";  // gray

      auto scaleColor = [&](float sc) -> const char*
      {
        if (!std::isfinite(sc) || sc < 0.0f) return C_WARN;
        if (sc > 1.05f) return C_GOOD;
        if (sc < 0.95f) return C_BAD;
        return C_WARN;
      };

      std::cout << "\n" << C_BAR << "====================================================================" << RST << "\n";
      std::cout << BOLD << C_INFO << "[JES PROBE] One-event before/after JetCalib (filtered)" << RST << "\n";

      std::cout << "  " << BOLD << C_RAW << "RAW" << RST
                << "   node: " << C_RAW << m_rawNode << RST << "  " << DIM << "(n=" << nRaw << ")" << RST << "\n";
      std::cout << "  " << BOLD << C_CAL << "CAL" << RST
                << "   node: " << C_CAL << m_calibNode << RST << "  " << DIM << "(n=" << nCalib << ")" << RST << "\n";

      std::cout << "  " << BOLD << "zvtx(MBD):" << RST << " " << C_INFO
                << std::fixed << std::setprecision(2) << zvtx << " cm" << RST
                << "  " << DIM << "(|z|<30 required)" << RST << "\n";

      std::cout << "  " << BOLD << "pt_raw >= " << RST << C_WARN
                << std::fixed << std::setprecision(1) << ptMinPrint << " GeV" << RST
                << "  " << DIM << "(printing " << nPass << " jets; maxRows=" << m_maxJets << ")" << RST << "\n";

      std::cout << C_BAR << "--------------------------------------------------------------------" << RST << "\n";

      // -------- fixed-width table layout (headers match rows) ----------
      const int W_I   = 3;
      const int W_PT  = 8;
      const int W_ETA = 8;
      const int W_PHI = 8;
      const int W_SC  = 6;

      // Header row (use SAME widths as data)
      std::cout
        << "  " << BOLD << std::setw(W_I) << "i" << RST << " " << C_BAR << "|" << RST << " "
        << BOLD << C_RAW << std::setw(W_PT)  << "pt_raw"  << RST << " "
        << BOLD << C_RAW << std::setw(W_ETA) << "eta_raw" << RST << " "
        << BOLD << C_RAW << std::setw(W_PHI) << "phi_raw" << RST << "  "
        << C_BAR << "=>" << RST << "  "
        << BOLD << C_CAL << std::setw(W_PT)  << "pt_cal"  << RST << " "
        << BOLD << C_CAL << std::setw(W_ETA) << "eta_cal" << RST << " "
        << BOLD << C_CAL << std::setw(W_PHI) << "phi_cal" << RST << "  "
        << C_BAR << "||" << RST << " "
        << BOLD << std::setw(W_SC) << "scale" << RST
        << "\n";

      // Separator line that matches the table width exactly
      auto rep = [](int n, char c) { return std::string(std::max(0, n), c); };
      const int tableWidth =
          2  // leading two spaces
        + W_I + 1 + 1 + 1  // i + space + '|' + space
        + W_PT + 1 + W_ETA + 1 + W_PHI  // raw cols with spaces
        + 2 + 2 + 2  // two spaces + "=>" + two spaces (approx; we print literal)
        + W_PT + 1 + W_ETA + 1 + W_PHI  // cal cols
        + 2 + 2 + 1 + W_SC; // two spaces + "||" + space + scale

      std::cout << C_BAR << rep(tableWidth, '-') << RST << "\n";

      int nPrinted = 0;
      for (int i = 0; i < nScan; ++i)
      {
        const Jet* jr = raw->get_jet(i);
        const Jet* jc = calib->get_jet(i);
        if (!jr || !jc) continue;

        const float pr = jr->get_pt();
        if (pr < ptMinPrint) continue;

        const float pc = jc->get_pt();
        const float sc = (pr > 0.0f) ? (pc / pr) : -1.0f;

        std::cout
          << "  " << BOLD << std::setw(W_I) << i << RST << " " << C_BAR << "|" << RST << " "
          << C_RAW
          << std::setw(W_PT)  << std::fixed << std::setprecision(3) << pr << " "
          << std::setw(W_ETA) << std::setprecision(3) << jr->get_eta() << " "
          << std::setw(W_PHI) << std::setprecision(3) << jr->get_phi()
          << RST << "  "
          << C_BAR << "=>" << RST << "  "
          << C_CAL
          << std::setw(W_PT)  << std::setprecision(3) << pc << " "
          << std::setw(W_ETA) << std::setprecision(3) << jc->get_eta() << " "
          << std::setw(W_PHI) << std::setprecision(3) << jc->get_phi()
          << RST << "  "
          << C_BAR << "||" << RST << " "
          << scaleColor(sc) << BOLD << std::setw(W_SC) << std::setprecision(3) << sc << RST
          << "\n";

        ++nPrinted;
        if (nPrinted >= m_maxJets) break;
      }

      std::cout << C_BAR << "====================================================================" << RST << "\n\n";

      m_fired = true;
      return Fun4AllReturnCodes::EVENT_OK;
  }

 private:
  std::string m_rawNode;
  std::string m_calibNode;
  int m_maxJets = 12;
  bool m_fired = false;
};



//======================================================================
//  The actual steering macro
//======================================================================
void Fun4All_recoilJets_unified_impl(const int   nEvents   =  0,
                                     const char* listFile  = "input_files.list",
                                     const char* outRoot   = "TrigPlot.root",
                                     const bool  verbose   = false)
{
    //--------------------------------------------------------------------
    // 0.  Banner & basic environment sanity
    //--------------------------------------------------------------------
    if (verbose) {
        std::cout << "\n>>> Fun4All_recoilJets – ana.495 driver <<<\n"
        << "    Input list : " << listFile  << '\n'
        << "    Output file: " << outRoot   << '\n'
        << "    nEvents    : " << nEvents   << (nEvents==0? " (all)\n":"\n");
    }
    
    Fun4AllServer* se = Fun4AllServer::instance();
    if (!se) detail::bail("unable to obtain Fun4AllServer instance!");
    
    
    auto env_lower = [](const char* key, const std::string& def = std::string{}) -> std::string
    {
        const char* raw = std::getenv(key);
        std::string s = raw ? detail::trim(std::string(raw)) : def;
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
        return s;
    };
    
    //--------------------------------------------------------------------
    // 1.  Parse the file list & determine run / segment
    //--------------------------------------------------------------------
    std::ifstream list(listFile);
    if (!list.is_open())
        detail::bail("cannot open input list \"" + std::string(listFile) + "\"");
    
    // ---------------------------------------------------------------
    // Parse list file
    //   DATA:  1 column  -> calo DST
    //   SIM :  5 columns -> calo + G4Hits + (truth jets) + global + mbd_epd
    // ---------------------------------------------------------------
    std::vector<std::string> filesCalo;
    std::vector<std::string> filesG4;
    std::vector<std::string> filesJets;
    std::vector<std::string> filesGlobal;
    std::vector<std::string> filesMbd;
    
    for (std::string line; std::getline(list, line); )
    {
        line = detail::trim(line);
        if (line.empty()) continue;
        if (!line.empty() && line[0] == '#') continue;
        
        // Columns:
        //   DATA : <DST_CALO_CLUSTER>
        //   SIM  : <DST_CALO_CLUSTER> <G4Hits> <DST_JETS> <DST_GLOBAL> <DST_MBD_EPD>
        std::istringstream iss(line);
        std::string fCalo, fG4, fJets, fGlobal, fMbd;
        iss >> fCalo >> fG4 >> fJets >> fGlobal >> fMbd;
        
        if (fCalo.empty()) continue;
        
        filesCalo.emplace_back(fCalo);
        
        // Keep vectors key-aligned (same length as filesCalo):
        // empty string == "not provided on this line"
        filesG4.emplace_back((!fG4.empty() && fG4 != "NONE") ? fG4 : std::string{});
        filesJets.emplace_back((!fJets.empty() && fJets != "NONE") ? fJets : std::string{});
        filesGlobal.emplace_back((!fGlobal.empty() && fGlobal != "NONE") ? fGlobal : std::string{});
        filesMbd.emplace_back((!fMbd.empty() && fMbd != "NONE") ? fMbd : std::string{});
    }
    
    if (filesCalo.empty())
        detail::bail("input list \"" + std::string(listFile) + "\" is empty");
    
    const std::string& firstFile = filesCalo.front(); // used for GetRunSegment
    
    // Dataset / input-mode detection
    //  - RJ_DATASET drives analysis mode: isPP | isPPrun25 | isAuAu | isSim | isSimEmbedded
    //  - RJ_CALO_INPUT_MODE optionally overrides input provenance: jetcalo | calofitting | simdst
    bool isSim = false;
    bool isSimEmbedded = false;
    bool isPPrun25 = false;
    bool isAuAuRequested = false;
    std::string datasetToken = env_lower("RJ_DATASET", "ispp");
    if (datasetToken == "issimembedded" || datasetToken == "simembedded")
    {
        isSim = true;
        isSimEmbedded = true;
        isAuAuRequested = true;
    }
    else if (datasetToken == "issim" || datasetToken == "sim"
             || datasetToken == "issimjet5" || datasetToken == "simjet5"
             || datasetToken == "issimmb" || datasetToken == "simmb")
    {
        isSim = true;
    }
    else if (datasetToken == "ispprun25" || datasetToken == "pprun25" || datasetToken == "pp25")
    {
        isPPrun25 = true;
    }
    else if (datasetToken == "isauau" || datasetToken == "auau" || datasetToken == "aa")
    {
        isAuAuRequested = true;
    }
    else
    {
        datasetToken = "ispp";
    }
    
    if (!isSim)
    {
        if (const char* f = std::getenv("RJ_IS_SIM"))
            isSim = (std::atoi(f) != 0);
        if (isSim)
        {
            datasetToken = "issim";
            isPPrun25 = false;
            isAuAuRequested = false;
        }
    }
    
    if (verbose)
    {
        std::cout << "[FLOW] dataset token resolution:"
        << " RJ_DATASET=" << datasetToken
        << " | isSim=" << (isSim ? "true" : "false")
        << " | isSimEmbedded=" << (isSimEmbedded ? "true" : "false")
        << " | isPPrun25=" << (isPPrun25 ? "true" : "false")
        << " | isAuAuRequested=" << (isAuAuRequested ? "true" : "false")
        << std::endl;
    }
    
    auto detect_input_mode = [&](const std::string& firstFileLower) -> std::string
    {
        std::string forced = env_lower("RJ_CALO_INPUT_MODE");
        if (forced == "jetcalo" || forced == "calofitting" || forced == "simdst") return forced;
        if (isSim && !isSimEmbedded) return "simdst";
        if (isSimEmbedded)
        {
            if (firstFileLower.find("dst_calo_") != std::string::npos) return "jetcalo";
            if (firstFileLower.find("calofitting") != std::string::npos) return "calofitting";
            return "jetcalo";
        }
        if (firstFileLower.find("calofitting") != std::string::npos) return "calofitting";
        if (firstFileLower.find("jetcalo") != std::string::npos) return "jetcalo";
        if (isPPrun25 || isAuAuRequested) return "calofitting";
        return "jetcalo";
    };
    
    std::string firstFileLower = firstFile;
    std::transform(firstFileLower.begin(), firstFileLower.end(), firstFileLower.begin(), [](unsigned char c){ return std::tolower(c); });
    std::string caloInputMode = detect_input_mode(firstFileLower);
    
    auto runSegPair = Fun4AllUtils::GetRunSegment(firstFile);
    int run = runSegPair.first;
    int seg = runSegPair.second;
    
    // MINIMAL / SCOPED FIX:
    // Fun4AllUtils::GetRunSegment can mis-read embedded Jet20 filenames like
    //   DST_CALO_Jet20-...-data-00054404-00001-00056.root
    // and return the embedded segment token as the "run".
    // For embedded/data calibration we need the real data run (>1000), so
    // repair that here whenever the parsed run is not a valid DATA timestamp.
    const std::size_t slashPos = firstFile.find_last_of("/\\");
    const std::string baseName = (slashPos == std::string::npos) ? firstFile : firstFile.substr(slashPos + 1);
    
    const bool useEmbeddedInclusiveJet20RunParse =
    isSimEmbedded &&
    baseName.compare(0, std::string("DST_CALO_Jet20-").size(), "DST_CALO_Jet20-") == 0 &&
    baseName.find("-data-") != std::string::npos;
    
    if (useEmbeddedInclusiveJet20RunParse && run <= 1000)
    {
        const std::string dataTag = "-data-";
        const std::size_t dataPos = baseName.find(dataTag);
        const std::string tail = baseName.substr(dataPos + dataTag.size());
        const std::size_t dash1 = tail.find('-');
        
        if (dash1 != std::string::npos)
        {
            const std::string runTok = tail.substr(0, dash1);
            const std::string rest = tail.substr(dash1 + 1);
            const std::size_t segEnd = rest.find_first_of("-.");
            const std::string segTok = rest.substr(0, segEnd);
            
            try
            {
                const int parsedRun = std::stoi(runTok);
                const int parsedSeg = std::stoi(segTok);
                if (parsedRun > 0)
                {
                    run = parsedRun;
                    seg = (parsedSeg >= 0) ? parsedSeg : 0;
                }
            }
            catch (...) { /* keep existing fallback logic below */ }
        }
    }
    
    if (run <= 0)
    {
        if (isSim)
        {
            // If the filename doesn't encode a run number, use a safe CDB timestamp.
            run = 1;
            seg = 0;
            
            // Optional override if you ever want it:
            //   export RJ_CDB_TIMESTAMP=XXXX
            if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
            {
                try
                {
                    const int t = std::stoi(ts);
                    if (t > 0) run = t;
                }
                catch (...) { /* keep default */ }
            }
            
            if (verbose)
                std::cout << "[INFO] Simulation mode: run/seg not in filename → using TIMESTAMP=" << run << "\n";
        }
        else
        {
            detail::bail("failed to extract run number from first file: " + firstFile);
        }
    }
    
    if (verbose)
        std::cout << "[INFO] Run=" << run << "  Seg=" << seg
        << "  (" << filesCalo.size() << " files)\n";
    
    // --------------------------------------------------------------------
    // Global verbosity control (RJ_VERBOSITY from env; defaults to 10;
    // Condor detection → 0). Also silences std::cout/cerr globally when 0.
    // --------------------------------------------------------------------
    int vlevel = 10;
    if (const char* venv = std::getenv("RJ_VERBOSITY"))
    {
        try { vlevel = std::stoi(venv); } catch (...) { vlevel = 10; }
    }
    else
    {
        // Detect Condor jobs defensively; default to quiet in batch.
        if (std::getenv("_CONDOR_SCRATCH_DIR") || std::getenv("_CONDOR_JOB_AD"))
            vlevel = 0;
    }
    
    // RAII silence for global std::cout/cerr if vlevel==0
    struct ScopedSilence {
        std::ofstream   sink;
        std::streambuf* cout_save = nullptr;
        std::streambuf* cerr_save = nullptr;
        bool active = false;
        void enable() {
            if (active) return;
            sink.open("/dev/null");
            cout_save = std::cout.rdbuf(sink.rdbuf());
            cerr_save = std::cerr.rdbuf(sink.rdbuf());
            active = true;
        }
        ~ScopedSilence() {
            if (active) {
                std::cout.rdbuf(cout_save);
                std::cerr.rdbuf(cerr_save);
            }
        }
    } _silence;
    
    if (vlevel == 0) _silence.enable();
    
    
    if (const char* fenv = std::getenv("RJ_F4A_VERBOSE"))
    {
        const int fv = std::atoi(fenv);
        if (fv > 0) se->Verbosity(fv);
    }
    
    // --------------------------------------------------------------------
    // YAML config (Phase-1): load once, validate lightly, and print summary
    // --------------------------------------------------------------------
    yamlcfg::Config cfg = yamlcfg::LoadConfig();
    if (const char* env = std::getenv("RJ_CLUSTER_UEPIPELINE"))
    {
        std::string s(env);
        if (s == "0" || s == "false" || s == "noSub") cfg.clusterUEpipeline = "noSub";
        else if (s == "1" || s == "true" || s == "variantA") cfg.clusterUEpipeline = "variantA";
        else if (s == "baseVariant") cfg.clusterUEpipeline = "baseVariant";
        else if (s == "variantB") cfg.clusterUEpipeline = "variantB";
        else cfg.clusterUEpipeline = s;  // pass through unknown for diagnostics
    }
    
    const std::vector<std::string> activeJetRKeys = yamlcfg::LoadJetRKeys(vlevel);
    
    std::vector<double> unfoldJetPtEdges;
    yamlcfg::ExpandUniformEdges(unfoldJetPtEdges,
                                cfg.unfold_jet_pt_start,
                                cfg.unfold_jet_pt_stop,
                                cfg.unfold_jet_pt_step);
    
    if (vlevel > 0)
    {
        std::cout << "\n[CFG] analysis_config.yaml\n"
        << "  path: " << cfg.yamlPath << "\n"
        << "  photon_eta_abs_max: " << cfg.photon_eta_abs_max << "\n"
        << "  jet_pt_min: " << cfg.jet_pt_min << "\n"
        << "  back_to_back_dphi_min_pi_fraction: " << cfg.back_to_back_dphi_min_pi_fraction
        << "  -> radians=" << (cfg.back_to_back_dphi_min_pi_fraction * M_PI) << "\n"
        << "  use_vz_cut: " << (cfg.use_vz_cut ? "true" : "false") << "\n"
        << "  vz_cut_cm: " << cfg.vz_cut_cm << "\n"
        << "  coneR: " << cfg.isoConeR << "\n"
        << "  matching: {pho_dr_max=" << cfg.pho_dr_max << ", jet_dr_max=" << cfg.jet_dr_max << "}\n"
        << "  isolation_wp: {aGeV=" << cfg.isoA << ", bPerGeV=" << cfg.isoB
        << ", sideGapGeV=" << cfg.isoGap << ", fixedGeV=" << cfg.isoFixed
        << ", coneR=" << cfg.isoConeR
        << ", towerMin=" << cfg.isoTowMin
        << ", isSlidingIso=" << (cfg.isSlidingIso ? "true" : "false") << "}\n"
        << "  isSlidingAndFixed: " << (cfg.isSlidingAndFixed ? "true" : "false") << "\n"
        << "  fixedGeV: " << cfg.isoFixed << "\n"
        << "  auau_cent_iso_wp: " << cfg.auauCentIsoWP.size() << " entries\n"
        << "  jes3_photon_pt_bins: [";
        for (std::size_t i = 0; i < cfg.jes3_photon_pt_bins.size(); ++i)
        {
            std::cout << cfg.jes3_photon_pt_bins[i] << (i + 1 < cfg.jes3_photon_pt_bins.size() ? ", " : "");
        }
        std::cout << "]\n  unfold_reco_photon_pt_bins: [";
        for (std::size_t i = 0; i < cfg.unfold_reco_photon_pt_bins.size(); ++i)
        {
            std::cout << cfg.unfold_reco_photon_pt_bins[i] << (i + 1 < cfg.unfold_reco_photon_pt_bins.size() ? ", " : "");
        }
        std::cout << "]\n  unfold_truth_photon_pt_bins: [";
        for (std::size_t i = 0; i < cfg.unfold_truth_photon_pt_bins.size(); ++i)
        {
            std::cout << cfg.unfold_truth_photon_pt_bins[i] << (i + 1 < cfg.unfold_truth_photon_pt_bins.size() ? ", " : "");
        }
        std::cout << "]\n  unfold_jet_pt_edges: start=" << cfg.unfold_jet_pt_start
        << " stop=" << cfg.unfold_jet_pt_stop
        << " step=" << cfg.unfold_jet_pt_step
        << " (nedges=" << unfoldJetPtEdges.size() << ")\n"
        << "  unfold_xj_bins: [";
        for (std::size_t i = 0; i < cfg.unfold_xj_bins.size(); ++i)
        {
            std::cout << cfg.unfold_xj_bins[i] << (i + 1 < cfg.unfold_xj_bins.size() ? ", " : "");
        }
        std::cout << "]\n"
        << "  clusterUEpipeline: " << cfg.clusterUEpipeline << "\n"
        << "  doPi0Analysis: " << (cfg.doPi0Analysis ? "true" : "false") << "\n"
        << "  event_display_tree: " << (cfg.event_display_tree ? "true" : "false") << "\n"
        << "  event_display_tree_max_per_bin: " << cfg.event_display_tree_max_per_bin << "\n\n";
        
        std::cout << "[CFG] caloInputMode: " << caloInputMode << "\n";
    }
    
    // --------------------------------------------------------------------
    // 2.  CDB + IO managers
    // --------------------------------------------------------------------
    recoConsts* rc = recoConsts::instance();
    // CDB_GLOBALTAG is REQUIRED for any CDBInterface::getUrl() call.
    // Allow override via env, otherwise default to your known-good tag.
    std::string gtag = "newcdbtag";
    if (const char* envGT = std::getenv("RJ_CDB_GLOBALTAG"))
    {
        std::string tmp = detail::trim(std::string(envGT));
        if (!tmp.empty()) gtag = tmp;
    }
    rc->set_StringFlag("CDB_GLOBALTAG", gtag);
    
    if (vlevel > 0)
        std::cout << "[INFO] CDB_GLOBALTAG=" << gtag << "\n";
    
    // TIMESTAMP:
    //  - DATA: use run number (must be > 1000 so Calo_Calib treats it as DATA)
    //  - SIM : use a known-good fixed timestamp (SIM does NOT run Process_Calo_Calib)
    unsigned long long cdbts = static_cast<unsigned long long>(run);
    
    if (!isSim || isSimEmbedded)
    {
        // DATA (and isSimEmbedded) must satisfy Calo_Calib's heuristic (TIMESTAMP > 1000 => data)
        // isSimEmbedded uses real AuAu calofitting DSTs, so it needs the actual run number.
        if (cdbts <= 1000ULL)
        {
            std::cerr << "[FATAL] DATA/embedded run number " << cdbts
            << " would make Calo_Calib treat this as SIM (TIMESTAMP<=1000).\n";
            throw std::runtime_error("Invalid DATA/embedded TIMESTAMP for Calo_Calib (must be > 1000).");
        }
    }
    else
    {
        cdbts = 47289ULL;  // keep your old working SIM timestamp
    }
    
    if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
    {
        char* end = nullptr;
        unsigned long long tmp = std::strtoull(ts, &end, 10);
        if (end != ts && tmp > 0ULL) cdbts = tmp;
    }
    
    rc->set_uint64Flag("TIMESTAMP", cdbts);
    
    if (vlevel > 0)
        std::cout << "[INFO] CDB TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP")
        << " (isSim=" << (isSim ? "true" : "false") << ")\n";
    
    CDBInterface::instance()->Verbosity(0);
    
    
    auto* flag = new FlagHandler();
    se->registerSubsystem(flag);
    
    // ------------------------------------------------------------------
    // Decide how to source truth jets in SIM:
    //
    //   RJ_TRUTH_JETS_MODE=AUTO  (default)
    //       - if list has 3rd column (DST_JETS): read jets directly from DST
    //       - otherwise: build truth jets from TRUTH particles (TruthJetInput)
    //
    //   RJ_TRUTH_JETS_MODE=DST
    //       - require 3rd column (DST_JETS) and read truth jets from DST
    //
    //   RJ_TRUTH_JETS_MODE=BUILD
    //       - ignore 3rd column and build truth jets from TRUTH particles
    //
    //   RJ_TRUTH_JETS_MODE=BOTH
    //       - read DST jets AND also build a second "from particles" truth-jet
    //         collection under a different node name for QA comparisons
    // ------------------------------------------------------------------
    bool useDSTTruthJets = false;
    bool buildTruthJetsFromParticles = false;
    bool buildTruthJetsAsAltNode = false;   // only true in BOTH mode
    
    auto all_nonempty = [](const std::vector<std::string>& v) -> bool
    {
        if (v.empty()) return false;
        for (const auto& s : v) if (s.empty()) return false;
        return true;
    };
    
    const bool listHasG4     = all_nonempty(filesG4);
    const bool listHasJets   = all_nonempty(filesJets);
    const bool listHasGlobal = all_nonempty(filesGlobal);
    const bool listHasMbd    = all_nonempty(filesMbd);
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] input contract:"
        << " caloFiles=" << filesCalo.size()
        << " | G4=" << (listHasG4 ? "present" : "missing")
        << " | DST_JETS=" << (listHasJets ? "present" : "missing")
        << " | DST_GLOBAL=" << (listHasGlobal ? "present" : "missing")
        << " | DST_MBD_EPD=" << (listHasMbd ? "present" : "missing")
        << " | caloInputMode=" << caloInputMode
        << std::endl;
    }
    
    // Normalize env token to lowercase
    
    const std::string truthMode = env_lower("RJ_TRUTH_JETS_MODE", "auto");
    if (truthMode == "dst")
    {
        useDSTTruthJets = true;
        buildTruthJetsFromParticles = false;
        buildTruthJetsAsAltNode = false;
    }
    else if (truthMode == "build")
    {
        useDSTTruthJets = false;
        buildTruthJetsFromParticles = true;
        buildTruthJetsAsAltNode = false;
    }
    else if (truthMode == "both")
    {
        useDSTTruthJets = true;
        buildTruthJetsFromParticles = true;
        buildTruthJetsAsAltNode = true;
    }
    else
    {
        // AUTO (or any unrecognized token): prefer DST truth jets if provided
        useDSTTruthJets = listHasJets;
        buildTruthJetsFromParticles = !listHasJets;
        buildTruthJetsAsAltNode = false;
    }
    
    // ------------------ Calo cluster DST (always) -------------------
    auto* inCalo = new Fun4AllDstInputManager("DSTcalofitting");
    for (const auto& f : filesCalo) inCalo->AddFile(f);
    se->registerInputManager(inCalo);
    
    if (isSim)
    {
        // For SIM we REQUIRE a reco-vertex stream from DST_GLOBAL.
        // DST_MBD_EPD is OPTIONAL: if absent (or 'NONE'), downstream code
        // will fall back to GlobalVertexMap for the reco vertex.
        // G4Hits is OPTIONAL: if missing, we skip truth-photon matching QA.
        if (!listHasGlobal)
        {
            std::ostringstream os;
            os << "isSim requires a reco-vertex DST_GLOBAL stream paired 1:1 with calo files.\n"
            << "Input list must include these columns per line:\n"
            << "  <DST_CALO_CLUSTER> <G4_or_NONE> <DST_JETS_or_NONE> <DST_GLOBAL> [<DST_MBD_EPD_or_NONE>]\n"
            << "But your list is missing DST_GLOBAL on at least one line.";
            detail::bail(os.str());
        }
        
        // G4 is OPTIONAL (photonjet productions may not provide it).
        // If you want to force it: export RJ_REQUIRE_G4=1
        bool requireG4 = false;
        if (const char* env = std::getenv("RJ_REQUIRE_G4")) requireG4 = (std::atoi(env) != 0);
        
        if (listHasG4)
        {
            auto* inG4 = isSimEmbedded
            ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_G4HITS_IN"))
            : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_G4HITS_IN"));
            for (const auto& f : filesG4) inG4->AddFile(f);
            se->registerInputManager(inG4);
        }
        else
        {
            if (requireG4)
            {
                detail::bail("RJ_REQUIRE_G4=1 but no G4Hits stream was provided in the input list.");
            }
            else
            {
                std::cout << "\033[31m[WARN] isSim: no G4Hits stream provided (or it's 'NONE'). "
                "Continuing; truth-photon matching QA will be skipped.\033[0m\n";
            }
        }
        
        auto* inGlobal = isSimEmbedded
        ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_GLOBAL_IN"))
        : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_GLOBAL_IN"));
        for (const auto& f : filesGlobal) inGlobal->AddFile(f);
        se->registerInputManager(inGlobal);
        
        if (listHasMbd)
        {
            auto* inMbd = isSimEmbedded
            ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_MBD_EPD_IN"))
            : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_MBD_EPD_IN"));
            for (const auto& f : filesMbd) inMbd->AddFile(f);
            se->registerInputManager(inMbd);
        }
        else
        {
            std::cout << "\033[31m[WARN] isSim: no DST_MBD_EPD stream provided (or it's 'NONE'). "
            "Continuing; reco vertex will use GlobalVertexMap fallback.\033[0m\n";
        }
        
        if (verbose)
            std::cout << "[INFO] isSim: registered input managers (Calo + Global"
            << (listHasMbd ? " + MBD_EPD" : " (no MBD_EPD)")
            << (listHasG4 ? " + G4" : " (no G4)") << ")\n";
    }
    
    
    // ------------------ Jets DST (SIM optional; for truth jets) -------
    if (isSim && useDSTTruthJets)
    {
        if (!listHasJets)
        {
            std::ostringstream os;
            os << "RJ_TRUTH_JETS_MODE=" << truthMode << " requires a list with at least 3 columns:\n"
            << "  <DST_CALO_CLUSTER> <G4Hits> <DST_JETS> [<DST_GLOBAL> <DST_MBD_EPD>]\n"
            << "Use your staged 5-column master list built from the matched lists.";
            detail::bail(os.str());
        }
        
        auto* inJets = isSimEmbedded
        ? static_cast<Fun4AllInputManager*>(new Fun4AllNoSyncDstInputManager("DST_JETS_IN"))
        : static_cast<Fun4AllInputManager*>(new Fun4AllDstInputManager("DST_JETS_IN"));
        for (const auto& f : filesJets) inJets->AddFile(f);
        se->registerInputManager(inJets);
        
        if (verbose)
            std::cout << "[INFO] isSim: registered DST_JETS input manager (truth jets from DST)\n";
    }
    
    if (verbose && isSim)
    {
        std::cout << "[INFO] Truth-jet mode (RJ_TRUTH_JETS_MODE=" << truthMode << "): "
        << (useDSTTruthJets ? "DST" : "")
        << ((useDSTTruthJets && buildTruthJetsFromParticles) ? "+" : "")
        << (buildTruthJetsFromParticles ? "BUILD" : "")
        << "\n";
    }
    
    // --------------------------------------------------------------------
    // 3.  Geometry + status + calibration + clustering
    // --------------------------------------------------------------------
    //
    // IMPORTANT:
    // RawClusterBuilderTemplate::process_event() ALWAYS requires "TOWERGEOM_CEMC".
    // But CaloGeomMapping with UseDetailedGeometry(true) publishes ONLY
    // "TOWERGEOM_CEMC_DETAILED" for CEMC.
    // So we ALWAYS create the legacy node (TOWERGEOM_CEMC), and optionally
    // also create the detailed node (TOWERGEOM_CEMC_DETAILED).
    //
    bool useDetailedCemcGeom = true;  // keep your current behavior by default
    if (const char* env = std::getenv("RJ_DETAILED_CEMC_GEOM"))
    {
        useDetailedCemcGeom = (std::atoi(env) != 0);
    }
    
    // Always publish the legacy/simple CEMC geometry node: TOWERGEOM_CEMC
    {
        auto* geomCemcLegacy = new CaloGeomMapping("Geom_CEMC");
        geomCemcLegacy->set_detector_name("CEMC");
        geomCemcLegacy->set_UseDetailedGeometry(false);
        se->registerSubsystem(geomCemcLegacy);
    }
    
    // Optionally publish the detailed CEMC geometry node: TOWERGEOM_CEMC_DETAILED
    if (useDetailedCemcGeom)
    {
        auto* geomCemcDetailed = new CaloGeomMapping("Geom_CEMC_DETAILED");
        geomCemcDetailed->set_detector_name("CEMC");
        geomCemcDetailed->set_UseDetailedGeometry(true);
        se->registerSubsystem(geomCemcDetailed);
    }
    
    // HCAL nodes (detailed not supported; CaloGeomMapping will fall back internally)
    for (const std::string& det : {"HCALIN","HCALOUT"})
    {
        auto* geom = new CaloGeomMapping(("Geom_" + det).c_str());
        geom->set_detector_name(det);
        geom->set_UseDetailedGeometry(true);
        se->registerSubsystem(geom);
    }
    
    
    // ------------------------------------------------------------------
    // Calo calibration + clustering contract
    //
    //   jetcalo     -> real-data lower-level tower input: run Process_Calo_Calib()
    //   calofitting -> real-data waveform-fit input: run Process_Calo_Calib()
    //   simdst      -> analysis DST already carries calibrated towers/clusters
    // ------------------------------------------------------------------
    if (isSim && !isSimEmbedded)
    {
        if (vlevel > 0)
        {
            std::cout << "[isSim] skipping Process_Calo_Calib() "
            << "(SIM DST already has TOWERINFO_CALIB and CLUSTERINFO_CEMC)\n";
        }
    }
    else if (caloInputMode == "calofitting" || caloInputMode == "jetcalo")
    {
        if (vlevel > 0)
        {
            if (isSimEmbedded)
            {
                if (caloInputMode == "calofitting")
                {
                    std::cout << "[isSimEmbedded] running Process_Calo_Calib() on CALOFITTING input\n";
                }
                else
                {
                    std::cout << "[isSimEmbedded] running Process_Calo_Calib() on DST_CALO / JETCALO input\n";
                }
            }
            else if (caloInputMode == "calofitting")
            {
                std::cout << "[DATA] running Process_Calo_Calib() on CALOFITTING input\n";
            }
            else
            {
                std::cout << "[DATA] running Process_Calo_Calib() on JETCALO input\n";
            }
            if (isAuAuRequested)
            {
                std::cout << "[DATA][AuAu] clusterUEpipeline may still apply native UE subtraction and reclusterization afterward\n";
            }
        }
        Process_Calo_Calib();
    }
    else
    {
        detail::bail("unsupported caloInputMode '" + caloInputMode + "'");
    }
    
    
    if (isSimEmbedded)
    {
        if (vlevel > 0)
        {
            std::cout << "[isSimEmbedded] skipping MbdReco (use embedded sample's existing MBD products)" << std::endl;
            std::cout << "[isSimEmbedded] skipping ZdcReco (not needed for embedded minimal path)" << std::endl;
            std::cout << "[isSimEmbedded] skipping GlobalVertexReco (use embedded sample's existing GlobalVertexMap)" << std::endl;
        }
    }
    else
    {
        if (vlevel > 0) std::cout << "Calibrating MBD" << std::endl;
        std::unique_ptr<MbdReco> mbdreco = std::make_unique<MbdReco>();
        se->registerSubsystem(mbdreco.release());
        
        if (!isSim && !isPPrun25)
        {
            if (vlevel > 0) std::cout << "Calibrating ZDC" << std::endl;
            auto* zdcreco = new ZdcReco();
            zdcreco->set_zdc1_cut(0.0);
            zdcreco->set_zdc2_cut(0.0);
            se->registerSubsystem(zdcreco);
        }
        else
        {
            if (vlevel > 0)
            {
                if (isPPrun25) std::cout << "[isPPrun25] skipping ZdcReco (CALOFITTING DST may not have TOWERS_ZDC)" << std::endl;
                else           std::cout << "[isSim] skipping ZdcReco (sim DST has no TOWERS_ZDC)" << std::endl;
            }
        }
        
        
        if (vlevel > 0) std::cout << "Retrieving Vtx Info" << std::endl;
        std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
        se->registerSubsystem(gvertex.release());
    }
    
    bool isAuAuData = false;
    bool isAuAuLike = false;
    if (isSimEmbedded)
    {
        isAuAuData = false;
        isAuAuLike = true;
    }
    else if (isSim)
    {
        isAuAuData = false;
        isAuAuLike = false;
    }
    else if (const char* env = std::getenv("RJ_DATASET"))
    {
        std::string sLower = detail::trim(std::string(env));
        std::transform(sLower.begin(), sLower.end(), sLower.begin(), [](unsigned char c){ return std::tolower(c); });
        if (sLower == "isauau" || sLower == "auau" || sLower == "aa")
            isAuAuData = true;
        else if (sLower == "ispp" || sLower == "pp" || sLower == "ispprun25" || sLower == "pprun25" || sLower == "pp25")
            isAuAuData = false;
        else
            isAuAuData = (run > 53864);
        isAuAuLike = isAuAuData;
    }
    else
    {
        isAuAuData = (run > 53864);
        isAuAuLike = isAuAuData;
    }
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] dataset semantics:"
        << " | isSim=" << (isSim ? "true" : "false")
        << " | isSimEmbedded=" << (isSimEmbedded ? "true" : "false")
        << " | isAuAuData=" << (isAuAuData ? "true" : "false")
        << " | isAuAuLike=" << (isAuAuLike ? "true" : "false")
        << " | isPPrun25=" << (isPPrun25 ? "true" : "false")
        << std::endl;
    }
    
    if (isAuAuLike && !isSimEmbedded)
    {
        if (vlevel > 0) std::cout << "building minbias classifier" << std::endl;
        auto* mb = new MinimumBiasClassifier();
        mb->Verbosity(0);
        mb->setOverwriteScale(
                              "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
        se->registerSubsystem(mb);
        
        if (vlevel > 0) std::cout << "building centrality classifier (Au+Au-like)" << std::endl;
        auto* cent = new CentralityReco();
        cent->Verbosity(0);
        cent->setOverwriteScale(
                                "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
        se->registerSubsystem(cent);
    }
    else
    {
        if (vlevel > 0)
        {
            if (isSimEmbedded) std::cout << "[isSimEmbedded] skipping MinimumBiasClassifier/CentralityReco (use embedded sample's existing centrality products)" << std::endl;
            else               std::cout << "[pp dataset] skipping CentralityReco" << std::endl;
        }
    }
    
    setenv("BEMCREC_CEMC_DISABLE_ASINH_POSITION", "0", 1);
    
    if (cfg.doPi0Analysis)
    {
        if (vlevel > 0)
        {
            std::cout << "[pi0] position-corrected-only mode: using CLUSTERINFO_CEMC only"
                      << " (no parallel CLUSTERINFO_CEMC_NOCORR branch will be built)" << std::endl;
        }
    }
    
    // ---------------------- Reco jets -----------------------------------------
    // Au+Au:
    //   - tower-level UE subtraction (RetowerCEMC + DTB + CASJ + DTB2 + SubtractTowers)
    //   - final jets are built from SUB1 tower containers
    //   - jets are written to: AntiKt_Tower_<rKey>_Sub1  (e.g. AntiKt_Tower_r04_Sub1)
    //
    // PP / isSim:
    //   - keep the existing pp-style JetReco + JetCalib chain (see 'else' below)
    // --------------------------------------------------------------------------
    std::string towerPrefixPCB = "TOWERINFO_CALIB";
    if (isAuAuLike)
    {
        if (const char* env = std::getenv("RJ_TOWERINFO_PREFIX"))
        {
            std::string s = detail::trim(std::string(env));
            if (!s.empty()) towerPrefixPCB = s;
        }
    }
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] reco/calibration branch:"
        << " | branch=" << (isAuAuLike ? "AuAu-like HI UE subtraction + SUB1 jets + JetCalib"
                            : "pp-style jets + JetCalib")
        << " | Process_Calo_Calib=" << (((!isSim || isSimEmbedded) && (caloInputMode == "calofitting" || caloInputMode == "jetcalo")) ? "ON" : "OFF")
        << " | clusterUEpipeline=" << cfg.clusterUEpipeline
        << " | towerPrefixPCB=" << towerPrefixPCB
        << " | truthJets=" << (useDSTTruthJets ? "DST" : "BUILD")
        << ((useDSTTruthJets && buildTruthJetsFromParticles) ? "+BUILD" : "")
        << std::endl;
    }
    
    if (isAuAuLike)
    {
        // Ensure JetBackground modules can attach nodes under DST/TOWER
        auto* ensure = new EnsureJetCalibNodes("EnsureJetCalibNodes_forHIJets");
        ensure->Verbosity(0);
        se->registerSubsystem(ensure);
        
        // Optional: control HI UE-subtraction verbosity independently
        int hiV = 0;
        if (const char* env = std::getenv("RJ_HIUE_VERBOSITY")) hiV = std::atoi(env);
        
        // Match Macro_HIJetReco.C semantics:
        //   0 = no flow
        //   1 = psi2 derived from calo
        //   2 = psi2 derived from HIJING
        //   3 = psi2 derived from sEPD
        int hiFlow = 0;
        if (const char* env = std::getenv("RJ_HI_DO_FLOW")) hiFlow = std::atoi(env);
        
        // TowerInfo node prefix for HI background chain.
        // calo-fitting Au+Au DSTs often publish per-detector nodes as:
        //   TOWERINFO_CEMC, TOWERINFO_HCALIN, TOWERINFO_HCALOUT
        // (NOT TOWERINFO_CALIB_*)
        std::string towerPrefix = "TOWERINFO_CALIB";
        if (const char* env = std::getenv("RJ_TOWERINFO_PREFIX"))
        {
            std::string s = detail::trim(std::string(env));
            if (!s.empty()) towerPrefix = s;
        }
        
        if (vlevel > 0)
            std::cout << "[HI] UE subtraction enabled: towerPrefix=" << towerPrefix
            << " do_flow=" << hiFlow
            << " (HIUE Verbosity=" << hiV << ")\n";
        
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
        if (hiFlow == 3)
        {
            auto* epreco = new EventPlaneReco();
            se->registerSubsystem(epreco);
        }
#endif
        
        // ------------------------------------------------------------------
        // 1) Retower CEMC (towerinfo)
        // ------------------------------------------------------------------
        auto* rcemc = new RetowerCEMC();
        rcemc->Verbosity(hiV);
        rcemc->set_towerinfo(true);
        rcemc->set_frac_cut(0.5);
        rcemc->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(rcemc);
        
        // ------------------------------------------------------------------
        // 2) Seed jets (RAW, R=0.2) for background estimation
        //     -> AntiKt_TowerInfo_HIRecoSeedsRaw_r02
        // ------------------------------------------------------------------
        {
            auto* seedReco = new JetReco("JetsReco_HIRecoSeedsRaw_r02");
            
            auto* incemc  = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, towerPrefix);
            auto* inihcal = new TowerJetInput(Jet::HCALIN_TOWERINFO,       towerPrefix);
            auto* inohcal = new TowerJetInput(Jet::HCALOUT_TOWERINFO,      towerPrefix);
            
            incemc->set_GlobalVertexType(GlobalVertex::MBD);
            inihcal->set_GlobalVertexType(GlobalVertex::MBD);
            inohcal->set_GlobalVertexType(GlobalVertex::MBD);
            
            seedReco->add_input(incemc);
            seedReco->add_input(inihcal);
            seedReco->add_input(inohcal);
            
            seedReco->add_algo(detail::fjAlgo(0.2f), "AntiKt_TowerInfo_HIRecoSeedsRaw_r02");
            seedReco->set_algo_node("ANTIKT");
            seedReco->set_input_node("TOWER");
            seedReco->Verbosity(hiV);
            se->registerSubsystem(seedReco);
        }
        
        // ------------------------------------------------------------------
        // 3) Background iteration 1 (seedType=0 -> RAW seeds)
        //     -> TowerInfoBackground_Sub1
        // ------------------------------------------------------------------
        auto* dtb = new DetermineTowerBackground();
        dtb->SetBackgroundOutputName("TowerInfoBackground_Sub1");
        dtb->SetFlow(hiFlow);
        dtb->SetSeedType(0);
        dtb->SetSeedJetD(3);
        dtb->Verbosity(hiV);
        dtb->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(dtb);
        
        // ------------------------------------------------------------------
        // 4) Copy + subtract jets using Sub1 background
        //     -> creates AntiKt_TowerInfo_HIRecoSeedsSub_r02
        // ------------------------------------------------------------------
        auto* casj = new CopyAndSubtractJets();
        casj->SetFlowModulation(hiFlow);
        casj->Verbosity(hiV);
        casj->set_towerinfo(true);
        casj->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(casj);
        
        // ------------------------------------------------------------------
        // 5) Background iteration 2 (seedType=1 -> SUB seeds)
        //     -> TowerInfoBackground_Sub2
        // ------------------------------------------------------------------
        auto* dtb2 = new DetermineTowerBackground();
        dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
        dtb2->SetFlow(hiFlow);
        dtb2->SetSeedType(1);
        dtb2->SetSeedJetPt(7);
        dtb2->Verbosity(hiV);
        dtb2->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(dtb2);
        
        // ------------------------------------------------------------------
        // 6) Subtract towers using Sub2 background
        //     -> writes *SUB1 tower containers used for final jet reco
        // ------------------------------------------------------------------
        auto* st = new SubtractTowers();
        st->SetFlowModulation(hiFlow);
        st->Verbosity(hiV);
        st->set_towerinfo(true);
        st->set_towerNodePrefix(towerPrefix);
        se->registerSubsystem(st);
        
        // ------------------------------------------------------------------
        // 7) Final jets from SUB1 tower containers (one node per R)
        //     -> RAW jets written to AntiKt_Tower_<rKey>_Sub1_RAW
        //     -> JetCalib output written to AntiKt_Tower_<rKey>_Sub1
        // ------------------------------------------------------------------
        int jetcalV = 0;
        if (const char* env = std::getenv("RJ_JETCALIB_VERBOSITY")) jetcalV = std::atoi(env);
        
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;
            
            const std::string calibNode = std::string("AntiKt_Tower_") + radKey + "_Sub1";
            const std::string rawNode   = calibNode + "_RAW";
            const std::string recoName  = std::string("JetsReco_AuAuSub_") + radKey;
            
            auto* jreco = new JetReco(recoName);
            
            auto* incemc  = new TowerJetInput(Jet::CEMC_TOWERINFO_SUB1,   towerPrefix);
            auto* inihcal = new TowerJetInput(Jet::HCALIN_TOWERINFO_SUB1, towerPrefix);
            auto* inohcal = new TowerJetInput(Jet::HCALOUT_TOWERINFO_SUB1, towerPrefix);
            
            incemc->set_GlobalVertexType(GlobalVertex::MBD);
            inihcal->set_GlobalVertexType(GlobalVertex::MBD);
            inohcal->set_GlobalVertexType(GlobalVertex::MBD);
            
            jreco->add_input(incemc);
            jreco->add_input(inihcal);
            jreco->add_input(inohcal);
            
            jreco->add_algo(detail::fjAlgo(R), rawNode);
            jreco->set_algo_node("ANTIKT");
            jreco->set_input_node("TOWER");
            
            int jetrecoV = 0;
            if (const char* env = std::getenv("RJ_JETRECO_VERBOSITY")) jetrecoV = std::atoi(env);
            jreco->Verbosity(jetrecoV);
            
            se->registerSubsystem(jreco);
            
            {
                auto* jcal = new JetCalib(std::string("JetCalib_AuAuSub_") + radKey);
                jcal->set_InputNode(rawNode);
                jcal->set_OutputNode(calibNode);
                jcal->set_JetRadius(R);
                jcal->set_ZvrtxNode("GlobalVertexMap");
                jcal->set_ApplyZvrtxDependentCalib(true);
                jcal->set_ApplyEtaDependentCalib(true);
                jcal->Verbosity(jetcalV);
                se->registerSubsystem(jcal);
                
                auto* probe = new JetCalibOneEventProbe(std::string("JetCalibOneEventProbe_AuAuSub_") + radKey,
                                                        rawNode,
                                                        calibNode,
                                                        /*maxJetsToPrint=*/12);
                probe->Verbosity(vlevel);
                se->registerSubsystem(probe);
                
                if (vlevel > 0)
                    std::cout << "[INFO] (AuAu) reco jets: built " << rawNode << " -> " << calibNode << " (R=" << R << ") from SUB1 towers with JetCalib\n";
            }
        }
    }
    else
    {
        // ---------------------- Reco jets + JES calibration (pp-style only) ----------------------
        //
        // IMPORTANT:
        // JetCalib::CreateNodeTree() requires a PHCompositeNode named "TOWER".
        // Many pp DSTs do NOT have it, so JetCalib aborts unless we create it.
        // We register EnsureJetCalibNodes once (only when we intend to run JetCalib).
        //
        // Apply pp JES calibration for ALL pp-style running (pp data AND isSim).
        // Keep Au+Au-like chains excluded.
        const bool doJetCalibAny = (!isAuAuLike);
        if (doJetCalibAny)
        {
            auto* ensure = new EnsureJetCalibNodes("EnsureJetCalibNodes_forJES");
            // keep this modest; set RJ_JETCALIB_NODE_VERBOSE=1 for prints
            int nodeV = 0;
            if (const char* env = std::getenv("RJ_JETCALIB_NODE_VERBOSE")) nodeV = std::atoi(env);
            ensure->Verbosity(nodeV);
            se->registerSubsystem(ensure);
            
            if (vlevel > 0)
                std::cout << "[INFO] JES: enabling JetCalib (pp-style: pp data + isSim) => ensuring DST/TOWER exists\n";
        }
        
        // Optional: control JetCalib verbosity independently
        int jetcalV = 0; // default: show InitRun/process_event messages
        if (const char* env = std::getenv("RJ_JETCALIB_VERBOSITY"))
        {
            jetcalV = std::atoi(env);
        }
        
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;
            
            // Canonical node name that RecoilJets reads (naming convention)
            const std::string calibNode = std::string("AntiKt_Tower_") + radKey;
            
            // Apply JES calibration for pp-like chains (pp data + pp-style SIM)
            // Run pp JES calibration in BOTH pp data and isSim (pp-style chains)
            const bool doJetCalib = (!isAuAuLike);
            
            // If calibrating: build RAW jets to a separate node to avoid name collision
            const std::string rawNode = doJetCalib ? (calibNode + "_RAW") : calibNode;
            
            // ---------------------- JetReco (build jets) ----------------------
            const std::string recoName = std::string("JetsReco_") + (isSim ? "Sim_" : "Data_") + radKey;
            
            auto* jreco = new JetReco(recoName);
            jreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO,    "TOWERINFO_CALIB"));
            jreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,  "TOWERINFO_CALIB"));
            jreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB"));
            
            jreco->add_algo(detail::fjAlgo(R), rawNode);
            jreco->set_algo_node("ANTIKT");
            jreco->set_input_node("TOWERINFO_CALIB");
            
            // If you want JetReco debug: export RJ_JETRECO_VERBOSITY=1
            int jetrecoV = 0;
            if (const char* env = std::getenv("RJ_JETRECO_VERBOSITY")) jetrecoV = std::atoi(env);
            jreco->Verbosity(jetrecoV);
            
            se->registerSubsystem(jreco);
            
            if (vlevel > 0)
            {
                std::cout << "[INFO] (" << (isSim ? "isSim" : "data") << ") reco jets: built "
                << rawNode << " (R=" << R << ") from TOWERINFO_CALIB\n";
            }
            
            // ---------------------- JetCalib (apply JES) ----------------------
            if (doJetCalib)
            {
                auto* jcal = new JetCalib(std::string("JetCalib_") + radKey);
                
                // JetCalib reads RAW jets and writes CALIB jets into the canonical node
                jcal->set_InputNode(rawNode);
                jcal->set_OutputNode(calibNode);
                
                jcal->set_JetRadius(R);
                jcal->set_ZvrtxNode("GlobalVertexMap");  // what JetCalib expects
                
                // Full pp JES: Zvrtx + eta dependent
                jcal->set_ApplyZvrtxDependentCalib(true);
                jcal->set_ApplyEtaDependentCalib(true);
                
                jcal->Verbosity(jetcalV);
                se->registerSubsystem(jcal);
                
                // ------------------------------------------------------------
                // Targeted JES probe: prints ONCE, ONLY if Verbosity() >= 20
                // ------------------------------------------------------------
                auto* probe = new JetCalibOneEventProbe(std::string("JetCalibOneEventProbe_") + radKey,
                                                        rawNode,
                                                        calibNode,
                                                        /*maxJetsToPrint=*/12);
                probe->Verbosity(vlevel);  // uses your RJ_VERBOSITY; probe prints only if >=20
                se->registerSubsystem(probe);
                
                if (vlevel > 0)
                {
                    std::cout << "[INFO] JES calib enabled: " << rawNode << " -> " << calibNode
                    << " (R=" << R << ", Zvrtx+eta dependent, JetCalibVerbosity=" << jetcalV << ")\n";
                }
            }
        }
    }
    
    // ---------------------- Truth jets -----------------------------------------
    // If useDSTTruthJets==true: they already exist from DST_JETS_IN.
    // If buildTruthJetsFromParticles==true: build them from TRUTH particles here.
    if (isSim && buildTruthJetsFromParticles)
    {
        if (vlevel > 0)
        {
            if (useDSTTruthJets)
                std::cout << "[INFO] (isSim) truth jets: ALSO building from TRUTH particles (QA 'both' mode)\n";
            else
                std::cout << "[INFO] (isSim) truth jets: building from TRUTH particles (no DST_JETS)\n";
        }
        
        for (const auto& radKey : activeJetRKeys)
        {
            int D = 0;
            try { D = std::stoi(radKey.substr(1)); } catch (...) { continue; }
            if (D <= 0) continue;
            const float R = 0.1f * D;
            
            auto* truthReco = new JetReco(std::string("TruthJetReco_FromParticles_") + radKey);
            auto* tji = new TruthJetInput(Jet::PARTICLE);
            tji->add_embedding_flag(1);  // pythia/herwig particles only
            truthReco->add_input(tji);
            
            // Node naming:
            //  - If we're NOT reading DST jets, keep canonical name "AntiKt_Truth_<rKey>"
            //  - If we ARE reading DST jets too (BOTH mode), avoid collisions
            const std::string truthNode = (useDSTTruthJets && buildTruthJetsAsAltNode)
            ? (std::string("AntiKt_TruthFromParticles_") + radKey)
            : (std::string("AntiKt_Truth_") + radKey);
            
            truthReco->add_algo(detail::fjAlgo(R), truthNode);
            truthReco->set_algo_node("ANTIKT");
            truthReco->set_input_node("TRUTH");
            
            // Same deal for truth-jet reco: keep JetReco quiet.
            truthReco->Verbosity(0);
            
            se->registerSubsystem(truthReco);
            
            if (vlevel > 0)
                std::cout << "[INFO] (isSim) truth jets: produced node " << truthNode
                << " (R=" << R << ")\n";
        }
    }
    else if (isSim && useDSTTruthJets && vlevel > 0)
    {
        std::cout << "[INFO] (isSim) truth jets: using nodes from DST_JETS (no TruthJetInput reco)\n";
    }
    
    
    // --------------------------------------------------------------------
    // 5.  Run-information helper (optional but handy)
    // --------------------------------------------------------------------
    if (!isSim)
    {
        auto* trigInfo = new TriggerRunInfoReco();
        trigInfo->Verbosity(vlevel);
        se->registerSubsystem(trigInfo);
    }
    else
    {
        if (vlevel > 0) std::cout << "[isSim] skipping TriggerRunInfoReco" << std::endl;
    }
    
    // Build photon clusters
    class NativeCEMCUESubtractor final : public SubsysReco
    {
    public:
        NativeCEMCUESubtractor(const std::string& name,
                               const std::string& inputNode,
                               const std::string& outputNode)
        : SubsysReco(name)
        , m_inputNode(inputNode)
        , m_outputNode(outputNode)
        {
        }
        
        void setDoSeedExclusion(bool v) { m_doSeedExclusion = v; }
        void setSeedJetNode(const std::string& n) { m_seedJetNode = n; }
        void setSeedMinPt(float v) { m_seedMinPt = v; }
        void setExclusionDR(float v) { m_exclusionDR = v; }
        
        int InitRun(PHCompositeNode* topNode) override
        {
            auto* src = findNode::getClass<TowerInfoContainer>(topNode, m_inputNode);
            if (!src)
            {
                std::cerr << Name() << ": missing input tower node '" << m_inputNode << "'" << std::endl;
                return Fun4AllReturnCodes::ABORTRUN;
            }
            
            PHNodeIterator iter(topNode);
            auto* cemcNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "CEMC"));
            if (!cemcNode)
            {
                std::cerr << Name() << ": missing CEMC composite node" << std::endl;
                return Fun4AllReturnCodes::ABORTRUN;
            }
            
            m_output = findNode::getClass<TowerInfoContainer>(topNode, m_outputNode);
            if (!m_output)
            {
                m_output = dynamic_cast<TowerInfoContainer*>(src->CloneMe());
                if (!m_output)
                {
                    std::cerr << Name() << ": failed to clone input tower container '" << m_inputNode << "'" << std::endl;
                    return Fun4AllReturnCodes::ABORTRUN;
                }
                
                auto* outNode = new PHIODataNode<PHObject>(m_output, m_outputNode, "PHObject");
                cemcNode->addNode(outNode);
            }
            
            return Fun4AllReturnCodes::EVENT_OK;
        }
        
        int process_event(PHCompositeNode* topNode) override
        {
            auto* src = findNode::getClass<TowerInfoContainer>(topNode, m_inputNode);
            auto* dst = findNode::getClass<TowerInfoContainer>(topNode, m_outputNode);
            
            if (!src || !dst)
            {
                std::cerr << Name()
                << ": missing required nodes (src=" << m_inputNode
                << ", dst=" << m_outputNode
                << ")"
                << std::endl;
                return Fun4AllReturnCodes::ABORTEVENT;
            }
            
            ++m_evt;
            
            // --- variantB: build seed-exclusion mask from coresoftware's refined seed jets ---
            std::vector<std::pair<float,float>> seedPositions;  // (eta, phi)
            RawTowerGeomContainer* geomCEMC_excl = nullptr;
            float exclusionDR2 = m_exclusionDR * m_exclusionDR;
            if (m_doSeedExclusion)
            {
                geomCEMC_excl = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
                auto* seeds = findNode::getClass<JetContainer>(topNode, m_seedJetNode);
                if (seeds && geomCEMC_excl)
                {
                    for (auto* jet : *seeds)
                    {
                        if (!jet) continue;
                        if (jet->get_pt() < m_seedMinPt) continue;
                        seedPositions.emplace_back(jet->get_eta(), jet->get_phi());
                    }
                }
                if (Verbosity() > 0)
                {
                    std::cout << "[" << Name() << "] evt=" << m_evt
                    << " seedExclusion: " << seedPositions.size()
                    << " seeds above " << m_seedMinPt << " GeV from " << m_seedJetNode
                    << " (DR=" << m_exclusionDR << ")" << std::endl;
                }
            }
            
            std::vector<double> stripSum;
            std::vector<unsigned int> stripCount;
            std::vector<double> stripMeanCache;
            
            double srcSumEAll = 0.0;
            double srcSumEGood = 0.0;
            unsigned int nSrcFiniteAll = 0;
            unsigned int nSrcFiniteGood = 0;
            unsigned int nSrcPositiveGood = 0;
            
            const unsigned int nchannels = src->size();
            for (unsigned int channel = 0; channel < nchannels; ++channel)
            {
                TowerInfo* srcTower = src->get_tower_at_channel(channel);
                if (!srcTower)
                {
                    continue;
                }
                
                const float energy = srcTower->get_energy();
                if (std::isfinite(energy))
                {
                    srcSumEAll += energy;
                    ++nSrcFiniteAll;
                }
                
                const unsigned int towerkey = src->encode_key(channel);
                const int ieta = src->getTowerEtaBin(towerkey);
                if (ieta < 0)
                {
                    continue;
                }
                
                if (static_cast<std::size_t>(ieta + 1) > stripSum.size())
                {
                    stripSum.resize(static_cast<std::size_t>(ieta + 1), 0.0);
                    stripCount.resize(static_cast<std::size_t>(ieta + 1), 0U);
                }
                
                if (!srcTower->get_isGood())
                {
                    continue;
                }
                
                if (!std::isfinite(energy))
                {
                    continue;
                }
                
                srcSumEGood += energy;
                ++nSrcFiniteGood;
                if (energy > 0.0f) ++nSrcPositiveGood;
                
                // variantB: skip towers within exclusion cone of any seed jet
                if (m_doSeedExclusion && !seedPositions.empty() && geomCEMC_excl)
                {
                    const int iphi = src->getTowerPhiBin(towerkey);
                    const RawTowerDefs::keytype gkey = RawTowerDefs::encode_towerid(
                                                                                    RawTowerDefs::CalorimeterId::CEMC, ieta, iphi);
                    RawTowerGeom* tg = geomCEMC_excl->get_tower_geometry(gkey);
                    if (tg)
                    {
                        const float tEta = tg->get_eta();
                        const float tPhi = tg->get_phi();
                        bool masked = false;
                        for (const auto& sp : seedPositions)
                        {
                            float deta = tEta - sp.first;
                            float dphi = tPhi - sp.second;
                            while (dphi >  M_PI) dphi -= 2.0f * M_PI;
                            while (dphi < -M_PI) dphi += 2.0f * M_PI;
                            if (deta * deta + dphi * dphi < exclusionDR2) { masked = true; break; }
                        }
                        if (masked) continue;
                    }
                }
                
                stripSum.at(static_cast<std::size_t>(ieta)) += energy;
                stripCount.at(static_cast<std::size_t>(ieta)) += 1U;
            }
            
            stripMeanCache.resize(stripSum.size(), 0.0);
            unsigned int nNonEmptyStrips = 0;
            double meanStripMean = 0.0;
            double maxStripMean = -std::numeric_limits<double>::infinity();
            int maxStripEta = -1;
            
            for (std::size_t ieta = 0; ieta < stripSum.size(); ++ieta)
            {
                const unsigned int nstrip = stripCount.at(ieta);
                if (nstrip == 0U) continue;
                
                const double stripMean = stripSum.at(ieta) / static_cast<double>(nstrip);
                stripMeanCache.at(ieta) = stripMean;
                meanStripMean += stripMean;
                ++nNonEmptyStrips;
                
                if (stripMean > maxStripMean)
                {
                    maxStripMean = stripMean;
                    maxStripEta = static_cast<int>(ieta);
                }
            }
            
            if (nNonEmptyStrips > 0)
            {
                meanStripMean /= static_cast<double>(nNonEmptyStrips);
            }
            else
            {
                maxStripMean = 0.0;
            }
            
            double dstSumEAll = 0.0;
            double dstSumEGood = 0.0;
            unsigned int nDstFiniteAll = 0;
            unsigned int nDstFiniteGood = 0;
            unsigned int nDstPositiveGood = 0;
            unsigned int nNegativeGood = 0;
            unsigned int nZeroedGood = 0;
            double mostNegative = 0.0;
            double maxSubtracted = -std::numeric_limits<double>::infinity();
            int maxSubtractedEta = -1;
            float maxSubtractedSrcE = 0.0f;
            float maxSubtractedDstE = 0.0f;
            
            for (unsigned int channel = 0; channel < nchannels; ++channel)
            {
                TowerInfo* srcTower = src->get_tower_at_channel(channel);
                TowerInfo* dstTower = dst->get_tower_at_channel(channel);
                if (!srcTower || !dstTower)
                {
                    continue;
                }
                
                const unsigned int towerkey = src->encode_key(channel);
                const int ieta = src->getTowerEtaBin(towerkey);
                
                float new_energy = 0.0f;
                float stripMean = 0.0f;
                const float src_energy = srcTower->get_energy();
                
                if (ieta >= 0 && static_cast<std::size_t>(ieta) < stripMeanCache.size() && srcTower->get_isGood())
                {
                    stripMean = static_cast<float>(stripMeanCache.at(static_cast<std::size_t>(ieta)));
                    if (std::isfinite(src_energy))
                    {
                        new_energy = src_energy - stripMean;
                    }
                }
                
                dstTower->set_time(srcTower->get_time());
                dstTower->set_energy(new_energy);
                
                if (std::isfinite(new_energy))
                {
                    dstSumEAll += new_energy;
                    ++nDstFiniteAll;
                }
                
                if (srcTower->get_isGood() && std::isfinite(new_energy))
                {
                    dstSumEGood += new_energy;
                    ++nDstFiniteGood;
                    if (new_energy > 0.0f) ++nDstPositiveGood;
                    if (new_energy < 0.0f)
                    {
                        ++nNegativeGood;
                        if (new_energy < mostNegative) mostNegative = new_energy;
                    }
                    if (std::fabs(new_energy) < 1e-6f) ++nZeroedGood;
                }
                
                if (srcTower->get_isGood() && std::isfinite(src_energy))
                {
                    const double subtracted = static_cast<double>(src_energy) - static_cast<double>(new_energy);
                    if (subtracted > maxSubtracted)
                    {
                        maxSubtracted = subtracted;
                        maxSubtractedEta = ieta;
                        maxSubtractedSrcE = src_energy;
                        maxSubtractedDstE = new_energy;
                    }
                }
            }
            
            if (Verbosity() > 0)
            {
                std::cout << "[" << Name() << "] evt=" << m_evt
                << " PHOSUB summary"
                << " | srcNode=" << m_inputNode
                << " dstNode=" << m_outputNode
                << " | nchannels=" << nchannels
                << " | strips(nonEmpty)=" << nNonEmptyStrips
                << " meanStrip=" << std::fixed << std::setprecision(4) << meanStripMean
                << " maxStrip=" << maxStripMean << "@ieta=" << maxStripEta
                << " | srcSumE(all/good)=" << std::setprecision(3) << srcSumEAll << "/" << srcSumEGood
                << " | dstSumE(all/good)=" << dstSumEAll << "/" << dstSumEGood
                << " | goodFinite(src/dst)=" << nSrcFiniteGood << "/" << nDstFiniteGood
                << " | goodE>0(src/dst)=" << nSrcPositiveGood << "/" << nDstPositiveGood
                << " | negGood=" << nNegativeGood
                << " zeroGood=" << nZeroedGood
                << " mostNegative=" << mostNegative
                << " | maxSubtracted=" << maxSubtracted
                << " (src=" << maxSubtractedSrcE
                << " -> dst=" << maxSubtractedDstE
                << ", ieta=" << maxSubtractedEta << ")"
                << std::endl;
            }
            
            if (Verbosity() > 1)
            {
                std::cout << "[" << Name() << "] evt=" << m_evt << " strip means:";
                int printed = 0;
                for (std::size_t ieta = 0; ieta < stripMeanCache.size(); ++ieta)
                {
                    if (stripCount.at(ieta) == 0U) continue;
                    std::cout << " (" << ieta
                    << ": n=" << stripCount.at(ieta)
                    << ", mean=" << std::fixed << std::setprecision(4) << stripMeanCache.at(ieta)
                    << ")";
                    ++printed;
                    if (printed >= 24)
                    {
                        std::cout << " ...";
                        break;
                    }
                }
                std::cout << std::endl;
            }
            
            return Fun4AllReturnCodes::EVENT_OK;
        }
        
    private:
        std::string m_inputNode;
        std::string m_outputNode;
        TowerInfoContainer* m_output{nullptr};
        int m_evt{0};
        bool m_doSeedExclusion{false};
        std::string m_seedJetNode{"AntiKt_TowerInfo_HIRecoSeedsSub_r02"};
        float m_seedMinPt{5.0f};
        float m_exclusionDR{0.4f};
    };
    
    class TowerInfoCanonicalRebaser final : public SubsysReco
    {
    public:
        TowerInfoCanonicalRebaser(const std::string& name,
                                  const std::string& sourceNode,
                                  const std::string& destNode)
        : SubsysReco(name)
        , m_sourceNode(sourceNode)
        , m_destNode(destNode)
        {
        }
        
        int process_event(PHCompositeNode* topNode) override
        {
            auto* src = findNode::getClass<TowerInfoContainer>(topNode, m_sourceNode);
            auto* dst = findNode::getClass<TowerInfoContainer>(topNode, m_destNode);
            
            if (!src || !dst)
            {
                std::cerr << Name()
                << ": missing source/destination tower nodes (src=" << m_sourceNode
                << ", dst=" << m_destNode << ")"
                << std::endl;
                return Fun4AllReturnCodes::ABORTEVENT;
            }
            
            const unsigned int nchannels = std::min(src->size(), dst->size());
            for (unsigned int channel = 0; channel < nchannels; ++channel)
            {
                TowerInfo* srcTower = src->get_tower_at_channel(channel);
                TowerInfo* dstTower = dst->get_tower_at_channel(channel);
                if (!srcTower || !dstTower)
                {
                    continue;
                }
                
                dstTower->set_time(srcTower->get_time());
                dstTower->set_energy(srcTower->get_energy());
            }
            
            return Fun4AllReturnCodes::EVENT_OK;
        }
        
    private:
        std::string m_sourceNode;
        std::string m_destNode;
    };
    
    std::string photonInputClusterNode = "CLUSTERINFO_CEMC";
    bool photonBuilderIsAuAu = false;
    
    if (cfg.clusterUEpipeline == "variantA" && isAuAuLike)
    {
        const std::string nativeCemcNode = towerPrefixPCB + "_CEMC_PHOSUB";
        
        int nativeUEV = 0;
        if (const char* env = std::getenv("RJ_HIUE_VERBOSITY")) nativeUEV = std::atoi(env);
        
        auto* nativeSub = new NativeCEMCUESubtractor("NativeCEMCUESubtractor",
                                                     towerPrefixPCB + "_CEMC",
                                                     nativeCemcNode);
        nativeSub->Verbosity(nativeUEV);
        se->registerSubsystem(nativeSub);
        
        if (nativeUEV > 0)
        {
            auto* auditBefore = new TowerAudit("TowerAudit_PHOSUB_before",
                                               towerPrefixPCB + "_CEMC",
                                               towerPrefixPCB + "_HCALIN_SUB1",
                                               towerPrefixPCB + "_HCALOUT_SUB1",
                                               10);
            auditBefore->Verbosity(nativeUEV);
            se->registerSubsystem(auditBefore);
            
            auto* auditAfter = new TowerAudit("TowerAudit_PHOSUB_after",
                                              nativeCemcNode,
                                              towerPrefixPCB + "_HCALIN_SUB1",
                                              towerPrefixPCB + "_HCALOUT_SUB1",
                                              10);
            auditAfter->Verbosity(nativeUEV);
            se->registerSubsystem(auditAfter);
        }
        
        // Recluster from UE-subtracted PHOSUB towers (ATLAS-like approach:
        // cluster on subtracted input so cluster energies, positions, and
        // tower membership all reflect the subtracted state)
        {
            auto* phoSubClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PHOSUB");
            phoSubClusterBuilder->Detector("CEMC");
            phoSubClusterBuilder->set_threshold_energy(0.070);
            const char* _calibroot = std::getenv("CALIBRATIONROOT");
            if (_calibroot && std::string(_calibroot).size())
            {
                std::string emc_prof_phosub = std::string(_calibroot) + "/EmcProfile/CEMCprof_Thresh30MeV.root";
                phoSubClusterBuilder->LoadProfile(emc_prof_phosub);
            }
            phoSubClusterBuilder->set_UseTowerInfo(1);
            phoSubClusterBuilder->set_UseAltZVertex(1);
            phoSubClusterBuilder->setInputTowerNodeName(nativeCemcNode);
            phoSubClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_PHOSUB");
            phoSubClusterBuilder->Verbosity(nativeUEV);
            se->registerSubsystem(phoSubClusterBuilder);
        }
        
        photonInputClusterNode = "CLUSTERINFO_CEMC_PHOSUB";
        
        if (vlevel > 0)
        {
            std::cout << "[clusterUEpipeline=variantA] enabled for AuAu-like running"
            << " | nativeCemcNode=" << nativeCemcNode
            << " | photonInputClusterNode=" << photonInputClusterNode
            << " | PhotonClusterBuilder will read PHOSUB clusters + PHOSUB/SUB1 tower nodes"
            << std::endl;
        }
    }
    else if (cfg.clusterUEpipeline == "variantB" && isAuAuLike)
    {
        const std::string nativeCemcNode = towerPrefixPCB + "_CEMC_PHOSUB";
        
        int nativeUEV = 0;
        if (const char* env = std::getenv("RJ_HIUE_VERBOSITY")) nativeUEV = std::atoi(env);
        
        auto* nativeSub = new NativeCEMCUESubtractor("NativeCEMCUESubtractor_B",
                                                     towerPrefixPCB + "_CEMC",
                                                     nativeCemcNode);
        nativeSub->setDoSeedExclusion(true);
        nativeSub->setSeedJetNode("AntiKt_TowerInfo_HIRecoSeedsSub_r02");
        nativeSub->setSeedMinPt(5.0f);
        nativeSub->setExclusionDR(0.4f);
        nativeSub->Verbosity(nativeUEV);
        se->registerSubsystem(nativeSub);
        
        // Recluster from seed-masked UE-subtracted PHOSUB towers
        {
            auto* phoSubClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate_PHOSUB");
            phoSubClusterBuilder->Detector("CEMC");
            phoSubClusterBuilder->set_threshold_energy(0.070);
            const char* _calibroot = std::getenv("CALIBRATIONROOT");
            if (_calibroot && std::string(_calibroot).size())
            {
                std::string emc_prof_phosub = std::string(_calibroot) + "/EmcProfile/CEMCprof_Thresh30MeV.root";
                phoSubClusterBuilder->LoadProfile(emc_prof_phosub);
            }
            phoSubClusterBuilder->set_UseTowerInfo(1);
            phoSubClusterBuilder->set_UseAltZVertex(1);
            phoSubClusterBuilder->setInputTowerNodeName(nativeCemcNode);
            phoSubClusterBuilder->setOutputClusterNodeName("CLUSTERINFO_CEMC_PHOSUB");
            phoSubClusterBuilder->Verbosity(nativeUEV);
            se->registerSubsystem(phoSubClusterBuilder);
        }
        
        photonInputClusterNode = "CLUSTERINFO_CEMC_PHOSUB";
        
        if (vlevel > 0)
        {
            std::cout << "[clusterUEpipeline=variantB] enabled for AuAu-like running"
            << " | nativeCemcNode=" << nativeCemcNode
            << " | seed exclusion: DR=0.4 from HIRecoSeedsSub_r02 (pt>5 GeV)"
            << " | photonInputClusterNode=" << photonInputClusterNode
            << std::endl;
        }
    }
    else if (cfg.clusterUEpipeline == "baseVariant" && isAuAuLike)
    {
        photonBuilderIsAuAu = true;
        
        if (vlevel > 0)
        {
            std::cout << "[clusterUEpipeline=baseVariant] enabled for AuAu-like running"
            << " | PCB m_is_auau=true → isolation uses RETOWER_SUB1/HCAL_SUB1 internally"
            << " | shower shapes use standard CEMC towers"
            << std::endl;
        }
    }
    
    auto* photonBuilder = new PhotonClusterBuilder("PhotonClusterBuilder");
    photonBuilder->set_input_cluster_node(photonInputClusterNode);
    photonBuilder->set_output_photon_node("PHOTONCLUSTER_CEMC");
    
    double minPhotonEt = 5.0;
    if (!cfg.jes3_photon_pt_bins.empty())
    {
        auto itMin = std::min_element(cfg.jes3_photon_pt_bins.begin(), cfg.jes3_photon_pt_bins.end());
        if (itMin != cfg.jes3_photon_pt_bins.end() && std::isfinite(*itMin))
        {
            minPhotonEt = *itMin;
        }
    }
    photonBuilder->set_ET_threshold(static_cast<float>(minPhotonEt));
    photonBuilder->set_iso_min_tower_energy(static_cast<float>(cfg.isoTowMin));
    
    photonBuilder->set_use_vz_cut(cfg.use_vz_cut);
    photonBuilder->set_vz_cut_cm(cfg.vz_cut_cm);
    
    photonBuilder->set_is_auau(photonBuilderIsAuAu);
    if ((cfg.clusterUEpipeline == "variantA" || cfg.clusterUEpipeline == "variantB") && isAuAuLike)
    {
        photonBuilder->set_emc_tower_node(towerPrefixPCB + "_CEMC_PHOSUB");
        photonBuilder->set_ihcal_tower_node(towerPrefixPCB + "_HCALIN_SUB1");
        photonBuilder->set_ohcal_tower_node(towerPrefixPCB + "_HCALOUT_SUB1");
    }
    else if (cfg.clusterUEpipeline == "baseVariant" && isAuAuLike)
    {
        photonBuilder->set_emc_tower_node(towerPrefixPCB + "_CEMC");
        photonBuilder->set_ihcal_tower_node(towerPrefixPCB + "_HCALIN");
        photonBuilder->set_ohcal_tower_node(towerPrefixPCB + "_HCALOUT");
        photonBuilder->set_tower_node_prefix(towerPrefixPCB);
    }
    else if (isAuAuLike)
    {
        // noSub: standard towers, no AuAu iso path
        photonBuilder->set_emc_tower_node(towerPrefixPCB + "_CEMC");
        photonBuilder->set_ihcal_tower_node(towerPrefixPCB + "_HCALIN");
        photonBuilder->set_ohcal_tower_node(towerPrefixPCB + "_HCALOUT");
    }
    
    if (vlevel > 0)
    {
        Dl_info pcbInfo{};
        if (dladdr((void*)&typeid(PhotonClusterBuilder), &pcbInfo) && pcbInfo.dli_fname)
            std::cout << "[DBG] PhotonClusterBuilder RTTI from: " << pcbInfo.dli_fname << "\n";
        else
            std::cout << "[DBG] PhotonClusterBuilder RTTI probe: dladdr failed\n";
        
        std::cout << "[DBG] PhotonClusterBuilder vzCut config: use="
        << (cfg.use_vz_cut ? "true" : "false")
        << " vz_cut_cm=" << cfg.vz_cut_cm
        << " | isAuAuLike=" << (isAuAuLike ? "true" : "false")
        << " | isSimEmbedded=" << (isSimEmbedded ? "true" : "false")
        << " | photonBuilderIsAuAu=" << (photonBuilderIsAuAu ? "true" : "false")
        << " | clusterUEpipeline=" << cfg.clusterUEpipeline
        << " | inputClusterNode=" << photonInputClusterNode << "\n";
    }
    
    
    photonBuilder->Verbosity(vlevel);
    se->registerSubsystem(photonBuilder);
    
    auto* recoilJets = new RecoilJets(outRoot);
    
    // ------------------------------------------------------------------
    // Apply YAML-driven knobs
    // ------------------------------------------------------------------
    recoilJets->setPhotonEtaAbsMax(cfg.photon_eta_abs_max);
    recoilJets->setMinJetPt(cfg.jet_pt_min);
    recoilJets->setMinBackToBack(cfg.back_to_back_dphi_min_pi_fraction * M_PI);
    
    recoilJets->setUseVzCut(cfg.use_vz_cut, cfg.vz_cut_cm);
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    recoilJets->setCentEdges(cfg.centrality_edges);
#endif
    recoilJets->setVertexReweighting(cfg.vertex_reweight_on,
                                     cfg.vertex_reweight_file,
                                     cfg.vertex_reweight_hist);
    recoilJets->setCentralityReweighting(cfg.centrality_reweight_on,
                                         cfg.centrality_reweight_file,
                                         cfg.centrality_reweight_hist);
    recoilJets->setActiveJetRKeys(activeJetRKeys);
    recoilJets->setIsolationWP(cfg.isoA, cfg.isoB, cfg.isoGap, cfg.isoConeR, cfg.isoTowMin, cfg.isoFixed);
    recoilJets->setIsSlidingIso(cfg.isSlidingIso);
    recoilJets->setTruthIsoMaxGeV(cfg.truthIsoGeV);
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
    if (!cfg.auauCentIsoWP.empty())
    {
        std::vector<RecoilJets::CentIsoWP> wps;
        wps.reserve(cfg.auauCentIsoWP.size());
        for (const auto& e : cfg.auauCentIsoWP)
            wps.push_back({e.aGeV, e.bPerGeV, e.sideGapGeV});
        recoilJets->setCentIsoWPs(wps);
    }
#endif
    
    recoilJets->setPhotonIDCuts(cfg.pre_e11e33_max,
                                cfg.pre_et1_min,
                                cfg.pre_et1_max,
                                cfg.pre_e32e35_min,
                                cfg.pre_e32e35_max,
                                cfg.pre_weta_max,
                                cfg.tight_w_lo,
                                cfg.tight_w_hi_intercept,
                                cfg.tight_w_hi_slope,
                                cfg.tight_e11e33_min,
                                cfg.tight_e11e33_max,
                                cfg.tight_et1_min,
                                cfg.tight_et1_max,
                                cfg.tight_e32e35_min,
                                cfg.tight_e32e35_max);
    
    recoilJets->setGammaPtBins(cfg.jes3_photon_pt_bins);
    recoilJets->setPhoMatchDRMax(cfg.pho_dr_max);
    recoilJets->setJetMatchDRMax(cfg.jet_dr_max);
    
    recoilJets->setUnfoldRecoPhotonPtBins(cfg.unfold_reco_photon_pt_bins);
    recoilJets->setUnfoldTruthPhotonPtBins(cfg.unfold_truth_photon_pt_bins);
    recoilJets->setUnfoldJetPtBins(unfoldJetPtEdges);
    recoilJets->setUnfoldXJBins(cfg.unfold_xj_bins);
    
    recoilJets->enablePi0Analysis(cfg.doPi0Analysis);
    recoilJets->setAnalysisConfigYAML(cfg.yamlText, "analysis_config.yaml");
    
    if (vlevel > 0)
    {
        std::cout << "[CFG] Applied to RecoilJets:"
        << " etaAbsMax=" << cfg.photon_eta_abs_max
        << " jetPtMin=" << cfg.jet_pt_min
        << " dphiMin(rad)=" << (cfg.back_to_back_dphi_min_pi_fraction * M_PI)
        << " useVzCut=" << (cfg.use_vz_cut ? "true" : "false")
        << " vzCut=" << cfg.vz_cut_cm
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
        << " centEdges_n=" << cfg.centrality_edges.size()
#endif
        << " phoDR=" << cfg.pho_dr_max
        << " jetDR=" << cfg.jet_dr_max
        << "\n";
        
        std::cout << "[CFG] isolation mode:";
        if (!cfg.isSlidingIso)
        {
            std::cout << " fixed (thrReco=fixedGeV=" << cfg.isoFixed << ")";
        }
        else
        {
#if defined(RJ_UNIFIED_ANALYSIS_AUAU)
            if (cfg.auauCentIsoWP.size() == 1)
            {
                std::cout << " sliding -> AuAu/embedded centrality-fit"
                << " (thrReco(cent)=" << cfg.auauCentIsoWP.front().aGeV
                << " + " << cfg.auauCentIsoWP.front().bPerGeV << " * cent)";
            }
            else if (!cfg.auauCentIsoWP.empty())
            {
                std::cout << " sliding -> per-centrality-bin WP list"
                << " (nCentWP=" << cfg.auauCentIsoWP.size() << ")";
            }
            else
#endif
            {
                std::cout << " sliding -> pT-dependent"
                << " (thrReco(pT)=" << cfg.isoA
                << " + " << cfg.isoB << " * pTgamma)";
            }
        }
        std::cout << " sideGap=" << cfg.isoGap << "\n";
    }
    
    recoilJets->enableEventDisplayDiagnostics(cfg.event_display_tree);
    recoilJets->setEventDisplayDiagnosticsMaxPerBin(cfg.event_display_tree_max_per_bin);
    
    if (vlevel > 0)
    {
        std::cout << "[CFG] EventDisplayTree: enable=" << (cfg.event_display_tree ? "true" : "false")
        << " max_per_bin=" << cfg.event_display_tree_max_per_bin << "\n";
    }
    
    recoilJets->Verbosity(vlevel);
    if (verbose) std::cout << "[INFO] RJ_VERBOSITY → " << vlevel << "\n";
    
    std::string dtype;
    if (isSim)
    {
        if (datasetToken == "issimembedded" || datasetToken == "simembedded")
            dtype = "isSimEmbedded";
        else if (datasetToken == "issimjet5" || datasetToken == "simjet5")
            dtype = "isSimJet5";
        else if (datasetToken == "issimmb" || datasetToken == "simmb")
            dtype = "isSimMB";
        else
            dtype = "isSim";
    }
    else if (isAuAuData) dtype = "isAuAu";
    else dtype = "isPP";
    if (isPPrun25) dtype = "isPP";
#if defined(RJ_UNIFIED_ANALYSIS_PP)
    if (dtype == "isAuAu" || dtype == "isSimEmbedded")
        detail::bail("Fun4All_recoilJets.C wrapper is pp analysis-module only; use Fun4All_recoilJets_AuAu.C for isAuAu/isSimEmbedded jobs.");
#endif
    
    if (verbose || vlevel > 0)
    {
        std::cout << "[FLOW] final RecoilJets mode:"
        << " dtype=" << dtype
        << " | datasetToken=" << datasetToken
        << " | photonInputClusterNode=" << photonInputClusterNode
        << " | photonBuilderIsAuAu=" << (photonBuilderIsAuAu ? "true" : "false")
        << std::endl;
    }
    recoilJets->setDataType(dtype);
    
    se->registerSubsystem(recoilJets);
    
    //--------------------------------------------------------------------
    // 6.  Run
    //--------------------------------------------------------------------
    try
    {
        if (vlevel > 0) std::cout << "[INFO] Starting event loop …" << std::endl;
        
        const bool stepEvents = ([]{
            const char* env = std::getenv("RJ_STEP_EVENTS");
            return (env && std::atoi(env) != 0);
        })();
        
        if (stepEvents && nEvents > 0)
        {
            for (int ievt = 0; ievt < nEvents; ++ievt)
            {
                std::cout << "[RUN] >>> event " << (ievt + 1) << "/" << nEvents << std::endl;
                const int rc = se->run(1);
                std::cout << "[RUN] <<< event " << (ievt + 1) << "/" << nEvents << "  rc=" << rc << std::endl;
                if (rc != 0) break;
            }
        }
        else
        {
            se->run(nEvents);
        }
        
        if (vlevel > 0) std::cout << "[INFO] Calling se->End() …" << std::endl;
        se->End();
        if (vlevel > 0) std::cout << "[INFO] Finished successfully." << std::endl;
    }
    catch (const std::exception& e)
    {
        detail::bail(std::string("exception in Fun4All: ") + e.what());
    }
}

#endif   // ROOT_VERSION guard

