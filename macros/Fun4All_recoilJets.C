//======================================================================
//  Fun4All_recoilJets.C
//  --------------------------------------------------------------------
#pragma once
#if defined(__CINT__) || defined(__CLING__)
  R__ADD_INCLUDE_PATH(/sphenix/u/patsfan753/thesisAnalysis/install/include)
#endif
#if defined(__CLING__)
  #pragma cling add_include_path("/sphenix/u/patsfan753/thesisAnalysis/install/include")
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

//–––– Standard Fun4All / sPHENIX ––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
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
#include <caloreco/PhotonClusterBuilder.h>
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
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cstdlib>     // getenv
#include <algorithm>   // std::transform
#include <cctype>      // std::tolower
#include <iomanip>     // std::setw, std::setprecision
#include <dlfcn.h>     // dlopen RTLD_NOLOAD
#include <TSystem.h>   // gSystem, GetBuildArch/Compiler info
#include <Calo_Calib.C>

// Load your local overrides FIRST so their symbols/dictionaries win
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libjetbase.so)


// Then load the rest of the environment stack
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libglobalvertex.so)
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

      double isoA = 1.08128;
      double isoB = 0.0299107;
      double isoGap = 1.0;
      double isoFixed = 2.0;
      double truthIsoGeV = 4.0;
      double isoConeR = 0.30;
      double isoTowMin = 0.0;
      bool   isSlidingIso = true;

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
void Fun4All_recoilJets(const int   nEvents   =  0,
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

    // Simulation detection:
    //  - primary: RJ_DATASET=isSim
    //  - fallback: RJ_IS_SIM=1 (wrapper sets this)
    //
    // CALOFITTING pp25 detection:
    //  - RJ_DATASET=isPPrun25 (or pp25/pprun25) enables the CALOFITTING calibration path
    bool isSim = false;
    bool isPPrun25 = false;
    if (const char* ds = std::getenv("RJ_DATASET"))
    {
        std::string s = detail::trim(std::string(ds));
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        if (s == "issim" || s == "sim") isSim = true;
        if (s == "ispprun25" || s == "pprun25" || s == "pp25") isPPrun25 = true;
      }
      if (!isSim)
      {
        if (const char* f = std::getenv("RJ_IS_SIM"))
          isSim = (std::atoi(f) != 0);
    }

  auto runSegPair = Fun4AllUtils::GetRunSegment(firstFile);
  int run = runSegPair.first;
  int seg = runSegPair.second;

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

    // --------------------------------------------------------------------
    // YAML config (Phase-1): load once, validate lightly, and print summary
    // --------------------------------------------------------------------
    yamlcfg::Config cfg = yamlcfg::LoadConfig();
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
        << "  matching: {pho_dr_max=" << cfg.pho_dr_max << ", jet_dr_max=" << cfg.jet_dr_max << "}\n"
        << "  isolation_wp: {aGeV=" << cfg.isoA << ", bPerGeV=" << cfg.isoB
        << ", sideGapGeV=" << cfg.isoGap << ", fixedGeV=" << cfg.isoFixed
        << ", truthIsoGeV=" << cfg.truthIsoGeV << ", coneR=" << cfg.isoConeR
        << ", towerMin=" << cfg.isoTowMin << "}\n"
        << "  isSlidingIso: " << (cfg.isSlidingIso ? "true" : "false") << "\n"
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
                  << "  event_display_tree: " << (cfg.event_display_tree ? "true" : "false") << "\n"
                  << "  event_display_tree_max_per_bin: " << cfg.event_display_tree_max_per_bin << "\n\n";
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

  if (!isSim)
  {
      // DATA must satisfy Calo_Calib's heuristic (TIMESTAMP > 1000 => data)
      if (cdbts <= 1000ULL)
      {
        std::cerr << "[FATAL] DATA run number " << cdbts
                  << " would make Calo_Calib treat this as SIM (TIMESTAMP<=1000).\n";
        throw std::runtime_error("Invalid DATA TIMESTAMP for Calo_Calib (must be > 1000).");
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

      // Normalize env token to lowercase
      auto env_lower = [](const char* key, const std::string& def) -> std::string
      {
        const char* raw = std::getenv(key);
        std::string s = raw ? detail::trim(std::string(raw)) : def;
        std::transform(s.begin(), s.end(), s.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return s;
      };

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
    auto* inCalo = new Fun4AllDstInputManager("DST_CALO_CLUSTER_IN");
    for (const auto& f : filesCalo) inCalo->AddFile(f);
    se->registerInputManager(inCalo);

    if (isSim)
    {
        // For SIM we REQUIRE reco-vertex streams (GLOBAL + MBD_EPD).
        // G4Hits is OPTIONAL: if missing, we skip truth-photon matching QA.
        if (!listHasGlobal || !listHasMbd)
        {
          std::ostringstream os;
          os << "isSim requires reco-vertex DSTs (DST_GLOBAL + DST_MBD_EPD) paired 1:1 with calo files.\n"
             << "Input list must include these columns per line:\n"
             << "  <DST_CALO_CLUSTER> <G4_or_NONE> <DST_JETS_or_NONE> <DST_GLOBAL> <DST_MBD_EPD>\n"
             << "But your list is missing DST_GLOBAL and/or DST_MBD_EPD on at least one line.";
          detail::bail(os.str());
        }

        // G4 is OPTIONAL (photonjet productions may not provide it).
        // If you want to force it: export RJ_REQUIRE_G4=1
        bool requireG4 = false;
        if (const char* env = std::getenv("RJ_REQUIRE_G4")) requireG4 = (std::atoi(env) != 0);

        if (listHasG4)
        {
          auto* inG4 = new Fun4AllDstInputManager("DST_G4HITS_IN");
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

        auto* inGlobal = new Fun4AllDstInputManager("DST_GLOBAL_IN");
        for (const auto& f : filesGlobal) inGlobal->AddFile(f);
        se->registerInputManager(inGlobal);

        auto* inMbd = new Fun4AllDstInputManager("DST_MBD_EPD_IN");
        for (const auto& f : filesMbd) inMbd->AddFile(f);
        se->registerInputManager(inMbd);

        if (verbose)
          std::cout << "[INFO] isSim: registered input managers (Calo + Global + MBD_EPD"
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

        auto* inJets = new Fun4AllDstInputManager("DST_JETS_IN");
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
    // Calo calibration + clustering
    //
    // DATA (isPP): run the standard Calo_Calib reconstruction chain.
    // DATA (isPPrun25): CALOFITTING DST → skip Process_Calo_Calib and instead:
    //   - CaloTowerCalib: TOWERS_* -> TOWERINFO_CALIB_*
    //   - RawClusterBuilderTemplate: ensure CLUSTERINFO_CEMC exists
    // SIM : skip (SIM analysis DST already has calibrated towers/clusters; re-running
    //             the reco chain is what triggers CaloTowerStatus hotmap/default-map failure).
    // ------------------------------------------------------------------
    if (!isSim && !isPPrun25)
    {
      if (vlevel > 0) std::cout << "[DATA] running Process_Calo_Calib()\n";
      Process_Calo_Calib();
    }
    else if (!isSim && isPPrun25)
    {
      if (vlevel > 0)
      {
        std::cout << "[isPPrun25] Skipping Process_Calo_Calib() "
                     "(CALOFITTING DST → re-calibrating towers & rebuilding clusters)\n";
        std::cout << "[isPPrun25] Running CaloTowerCalib: inputPrefix=TOWERS_ -> outputPrefix=TOWERINFO_CALIB_\n";
      }

      CaloTowerCalib* calibEMC = new CaloTowerCalib("CaloTowerCalib_CEMC_fromTOWERS");
      calibEMC->set_detector_type(CaloTowerDefs::CEMC);
      calibEMC->set_inputNodePrefix("TOWERS_");
      calibEMC->set_outputNodePrefix("TOWERINFO_CALIB_");
      calibEMC->set_doCalibOnly(true);
      se->registerSubsystem(calibEMC);

      CaloTowerCalib* calibIHCal = new CaloTowerCalib("CaloTowerCalib_HCALIN_fromTOWERS");
      calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
      calibIHCal->set_inputNodePrefix("TOWERS_");
      calibIHCal->set_outputNodePrefix("TOWERINFO_CALIB_");
      calibIHCal->set_doCalibOnly(true);
      se->registerSubsystem(calibIHCal);

      CaloTowerCalib* calibOHCal = new CaloTowerCalib("CaloTowerCalib_HCALOUT_fromTOWERS");
      calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
      calibOHCal->set_inputNodePrefix("TOWERS_");
      calibOHCal->set_outputNodePrefix("TOWERINFO_CALIB_");
      calibOHCal->set_doCalibOnly(true);
      se->registerSubsystem(calibOHCal);

      if (vlevel > 0) std::cout << "[isPPrun25] Building clusters: RawClusterBuilderTemplate -> CLUSTERINFO_CEMC\n";

      RawClusterBuilderTemplate* ClusterBuilder =
          new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
      ClusterBuilder->Detector("CEMC");
      ClusterBuilder->set_threshold_energy(0.070);  // match Process_Calo_Calib default

      const char* calibroot = std::getenv("CALIBRATIONROOT");
      if (calibroot && std::string(calibroot).size())
      {
        std::string emc_prof = std::string(calibroot) + "/EmcProfile/CEMCprof_Thresh30MeV.root";
        ClusterBuilder->LoadProfile(emc_prof);
        if (vlevel > 0) std::cout << "[isPPrun25] ClusterBuilder LoadProfile: " << emc_prof << "\n";
      }
      else
      {
        if (vlevel > 0) std::cout << "[isPPrun25][WARN] CALIBRATIONROOT not set; cluster profile not loaded\n";
      }

      ClusterBuilder->set_UseTowerInfo(1);
      ClusterBuilder->set_UseAltZVertex(1);
      se->registerSubsystem(ClusterBuilder);
    }
    else
    {
      if (vlevel > 0)
      {
        std::cout << "[isSim] skipping Process_Calo_Calib() "
                     "(SIM DST already has TOWERINFO_CALIB and CLUSTERINFO_CEMC)\n";
      }
    }


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

  // pp-only steering macro: Au+Au functionality is handled in Fun4All_recoilJets_AuAu.C
  // Keep behavior identical for pp data and isSim (pp-style chain).
  const bool isAuAuData = false;

  if (vlevel > 0)
      std::cout << "[pp macro] forcing pp-style running (no Au+Au centrality/minbias modules)\n";

  // ---------------------- Reco jets + JES calibration (pp-style only) ----------------------
  //
  // IMPORTANT:
  // JetCalib::CreateNodeTree() requires a PHCompositeNode named "TOWER".

  // Apply pp JES calibration for ALL pp-style running (pp data AND isSim).
  // Keep Au+Au excluded.
  const bool doJetCalibAny = (!isAuAuData);
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
      const bool doJetCalib = (!isAuAuData);

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
  auto* photonBuilder = new PhotonClusterBuilder("PhotonClusterBuilder");
  photonBuilder->set_input_cluster_node("CLUSTERINFO_CEMC");
  photonBuilder->set_output_photon_node("PHOTONCLUSTER_CEMC");
  photonBuilder->set_vz_cut(cfg.use_vz_cut, cfg.vz_cut_cm);
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
  recoilJets->setActiveJetRKeys(activeJetRKeys);
  recoilJets->setIsolationWP(cfg.isoA, cfg.isoB, cfg.isoGap, cfg.isoConeR, cfg.isoTowMin, cfg.isoFixed);
  recoilJets->setIsSlidingIso(cfg.isSlidingIso);
  recoilJets->setTruthIsoMaxGeV(cfg.truthIsoGeV);

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

  recoilJets->setAnalysisConfigYAML(cfg.yamlText, "analysis_config.yaml");

  if (vlevel > 0)
  {
      std::cout << "[CFG] Applied to RecoilJets:"
                << " etaAbsMax=" << cfg.photon_eta_abs_max
                << " jetPtMin=" << cfg.jet_pt_min
                << " dphiMin(rad)=" << (cfg.back_to_back_dphi_min_pi_fraction * M_PI)
                << " useVzCut=" << (cfg.use_vz_cut ? "true" : "false")
                << " vzCut=" << cfg.vz_cut_cm
                << " phoDR=" << cfg.pho_dr_max
                << " jetDR=" << cfg.jet_dr_max
                << "\n";
  }
      
    // EventDisplay diagnostics payload (EventDisplayTree written into the ROOT output)
    recoilJets->enableEventDisplayDiagnostics(cfg.event_display_tree);
    recoilJets->setEventDisplayDiagnosticsMaxPerBin(cfg.event_display_tree_max_per_bin);

    if (vlevel > 0)
    {
      std::cout << "[CFG] EventDisplayTree: enable=" << (cfg.event_display_tree ? "true" : "false")
                << " max_per_bin=" << cfg.event_display_tree_max_per_bin << "\n";
    }

    
  // RecoilJets inherits SubsysReco::Verbosity(int)
  recoilJets->Verbosity(vlevel);
  if (verbose) std::cout << "[INFO] RJ_VERBOSITY → " << vlevel << '\n';
  // Pick analysis type for the module (pp macro supports only isPP / isSim).
  std::string dtype = "isPP";
  if (const char* env = std::getenv("RJ_DATASET"))
  {
              std::string s = detail::trim(std::string(env));
              std::string sLower = s;
              std::transform(sLower.begin(), sLower.end(), sLower.begin(),
                             [](unsigned char c){ return std::tolower(c); });

              if (sLower == "issim" || sLower == "sim")      dtype = "isSim";
              else                                           dtype = "isPP";
  }

  if (verbose) std::cout << "[INFO] RJ_DATASET → " << dtype << '\n';
  recoilJets->setDataType(dtype);

  se->registerSubsystem(recoilJets);

  //--------------------------------------------------------------------
  // 6.  Run
  //--------------------------------------------------------------------
  try
  {
    if (verbose) std::cout << "[INFO] Starting event loop …\n";
    se->run(nEvents);
    se->End();
    if (verbose) std::cout << "[INFO] Finished successfully.\n";
  }
  catch (const std::exception& e)
  {
    detail::bail(std::string("exception in Fun4All: ") + e.what());
  }

//  //--------------------------------------------------------------------
//  // 7.  Clean exit
//  //--------------------------------------------------------------------
//  gSystem->Exit(0);
}

#endif   // ROOT_VERSION guard
