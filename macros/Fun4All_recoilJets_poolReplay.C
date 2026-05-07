#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TKey.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TObjString.h>
#include <TProfile.h>
#include <TProfile3D.h>
#include <TTree.h>
#include <TVector2.h>
#include <TMVA/RBDT.hxx>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace rjpool
{
  constexpr double PRE_E11E33_MAX = 0.98;
  constexpr double PRE_ET1_MIN = 0.60;
  constexpr double PRE_ET1_MAX = 1.00;
  constexpr double PRE_E32E35_MIN = 0.80;
  constexpr double PRE_E32E35_MAX = 1.00;
  constexpr double PRE_WETA_MAX = 0.60;

  constexpr double TIGHT_W_LO = 0.0;
  constexpr double TIGHT_E11E33_MIN = 0.40;
  constexpr double TIGHT_E11E33_MAX = 0.98;
  constexpr double TIGHT_ET1_MIN = 0.90;
  constexpr double TIGHT_ET1_MAX = 1.00;
  constexpr double TIGHT_E32E35_MIN = 0.92;
  constexpr double TIGHT_E32E35_MAX = 1.00;

  std::string trim(std::string s)
  {
    auto notSpace = [](unsigned char c) { return !std::isspace(c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
    s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());
    return s;
  }

  bool startsWithKey(const std::string& line, const std::string& key)
  {
    const std::string p = key + ":";
    return line.rfind(p, 0) == 0;
  }

  std::string afterColon(const std::string& line)
  {
    const std::size_t p = line.find(':');
    return p == std::string::npos ? std::string{} : trim(line.substr(p + 1));
  }

  bool parseBool(std::string s, bool& out)
  {
    s = trim(s);
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    if (s == "true" || s == "1" || s == "yes" || s == "on") { out = true; return true; }
    if (s == "false" || s == "0" || s == "no" || s == "off") { out = false; return true; }
    return false;
  }

  bool parseDouble(std::string s, double& out)
  {
    s = trim(s);
    const std::size_t hash = s.find('#');
    if (hash != std::string::npos) s = trim(s.substr(0, hash));
    if (!s.empty() && s.front() == '[')
    {
      s.erase(s.begin());
      const std::size_t comma = s.find(',');
      const std::size_t close = s.find(']');
      std::size_t end = std::string::npos;
      if (comma != std::string::npos && close != std::string::npos) end = std::min(comma, close);
      else if (comma != std::string::npos) end = comma;
      else if (close != std::string::npos) end = close;
      if (end != std::string::npos) s = s.substr(0, end);
    }
    s = trim(s);
    if (s.empty()) return false;
    char* endptr = nullptr;
    const double v = std::strtod(s.c_str(), &endptr);
    if (endptr == s.c_str()) return false;
    out = v;
    return true;
  }

  void parseDoubles(std::string s, std::vector<double>& out)
  {
    out.clear();
    const std::size_t l = s.find('[');
    const std::size_t r = s.find(']');
    if (l == std::string::npos || r == std::string::npos || r <= l) return;
    std::stringstream ss(s.substr(l + 1, r - l - 1));
    for (std::string tok; std::getline(ss, tok, ','); )
    {
      double v = 0.0;
      if (parseDouble(tok, v)) out.push_back(v);
    }
  }

  void parseInts(std::string s, std::vector<int>& out)
  {
    out.clear();
    std::vector<double> vals;
    parseDoubles(s, vals);
    for (double v : vals) out.push_back(static_cast<int>(std::lround(v)));
  }

  void parseStrings(std::string s, std::vector<std::string>& out)
  {
    out.clear();
    const std::size_t hash = s.find('#');
    if (hash != std::string::npos) s = s.substr(0, hash);
    const std::size_t l = s.find('[');
    const std::size_t r = s.rfind(']');
    if (l == std::string::npos || r == std::string::npos || r <= l) return;
    std::stringstream ss(s.substr(l + 1, r - l - 1));
    for (std::string tok; std::getline(ss, tok, ','); )
    {
      tok = trim(tok);
      if (tok.size() >= 2 &&
          ((tok.front() == '"' && tok.back() == '"') ||
           (tok.front() == '\'' && tok.back() == '\'')))
      {
        tok = tok.substr(1, tok.size() - 2);
      }
      tok = trim(tok);
      if (!tok.empty()) out.push_back(tok);
    }
  }

  bool hasFeature(const std::vector<std::string>& features, const std::string& wanted)
  {
    for (const std::string& f : features)
    {
      std::string k = f;
      std::transform(k.begin(), k.end(), k.begin(), [](unsigned char c) { return std::tolower(c); });
      if (k == wanted) return true;
    }
    return false;
  }

  std::vector<std::string> withCentralityFeature(std::vector<std::string> features)
  {
    if (!hasFeature(features, "centrality") && !hasFeature(features, "cent"))
      features.push_back("centrality");
    return features;
  }

  bool parseInlineMapDouble(const std::string& line, const std::string& key, double& out)
  {
    const std::size_t p = line.find(key + ":");
    if (p == std::string::npos) return false;
    std::size_t start = p + key.size() + 1;
    while (start < line.size() && std::isspace(static_cast<unsigned char>(line[start]))) ++start;
    std::size_t end = start;
    while (end < line.size() && line[end] != ',' && line[end] != '}') ++end;
    return parseDouble(line.substr(start, end - start), out);
  }

  std::string normalizePre(std::string s)
  {
    s = trim(s);
    std::string k = s;
    std::transform(k.begin(), k.end(), k.begin(), [](unsigned char c) { return std::tolower(c); });
    if (k.empty() || k == "reference") return "reference";
    if (k == "varianta" || k == "newppg12") return "newPPG12";
    if (k == "variantb" || k == "noprecriteria") return "noPreCriteria";
    if (k == "variantc" || k == "onlynpb") return "onlyNPB";
    if (k == "variantd" || k == "refplusnpb") return "refPlusNPB";
    if (k == "variante" || k == "auauonlynpb") return "auauOnlyNPB";
    return s;
  }

  std::string normalizeTight(std::string s)
  {
    s = trim(s);
    std::string k = s;
    std::transform(k.begin(), k.end(), k.begin(), [](unsigned char c) { return std::tolower(c); });
    if (k.empty() || k == "reference") return "reference";
    if (k == "varianta" || k == "newppg12") return "newPPG12";
    if (k == "variantb" || k == "auauembeddedbdt") return "auauEmbeddedBDT";
    if (k == "centindcontrol") return "centINDcontrol";
    if (k == "centasfeat" || k == "centasfeature") return "centAsFeat";
    if (k == "centdepbdts" || k == "centdepbdt") return "centDepBDTs";
    return s;
  }

  std::string normalizeNonTight(std::string s)
  {
    s = trim(s);
    std::string k = s;
    std::transform(k.begin(), k.end(), k.begin(), [](unsigned char c) { return std::tolower(c); });
    if (k.empty() || k == "reference") return "reference";
    if (k == "varianta" || k == "newppg12" || k == "bdtsideband") return "newPPG12";
    if (k == "variantb" || k == "auaubdtsideband") return "auauBDTSideband";
    if (k == "variantc" || k == "auaubdtcomplement") return "auauBDTComplement";
    if (k == "centindcontrol") return "centINDcontrol";
    if (k == "centasfeat" || k == "centasfeature") return "centAsFeat";
    if (k == "centdepbdts" || k == "centdepbdt") return "centDepBDTs";
    return s;
  }

  struct Config
  {
    bool useVzCut = true;
    double vzCutCm = 60.0;
    std::vector<double> ptBins = {15,17,19,21,23,26,35};
    std::vector<double> unfoldRecoPhotonPtBins = {10,15,17,19,21,23,26,35,40};
    std::vector<double> unfoldTruthPhotonPtBins = {5,10,15,17,19,21,23,26,35,40};
    std::vector<double> unfoldXJBins = {0.0,0.20,0.24,0.29,0.35,0.41,0.50,0.60,0.72,0.86,1.03,1.24,1.49,1.78,2.14,3.0};
    std::vector<double> jetRadii = {0.2, 0.3, 0.4};
    std::vector<int> centEdges = {0,20,50,80};
    double minJetPt = 10.0;
    double minBackToBack = 7.0 * M_PI / 8.0;
    double coneR = 0.4;
    bool slidingIso = false;
    double isoA = 1.08128;
    double isoB = 0.0299107;
    double isoGap = 1.0;
    double isoFixed = 2.0;
    double npbCut = 0.5;
    double auauNpbCut = 0.5;
    double tightBDTMinIntercept = 0.8333333333333334;
    double tightBDTMinSlope = -0.003333333333333336;
    double nonTightBDTMinIntercept = 0.7333333333333333;
    double nonTightBDTMinSlope = -0.01333333333333333;
    double nonTightBDTMaxIntercept = 0.6666666666666666;
    double nonTightBDTMaxSlope = 0.003333333333333336;
    double auauTightBDTMinIntercept = 0.8333333333333334;
    double auauTightBDTMinSlope = -0.003333333333333336;
    double auauTightBDTMax = 1.0;
    double auauNonTightBDTMinIntercept = 0.7333333333333333;
    double auauNonTightBDTMinSlope = -0.01333333333333333;
    double auauNonTightBDTMaxIntercept = 0.6666666666666666;
    double auauNonTightBDTMaxSlope = 0.003333333333333336;
    std::string auauNpbModelFile;
    std::vector<std::string> auauNpbFeatures;
    std::string auauTightBDTModelFile;
    std::string auauTightBDTCentINDControlModelFile;
    std::string auauTightBDTCentAsFeatModelFile;
    std::vector<std::string> auauTightBDTCentDepModelFiles;
    std::vector<std::string> auauTightBDTFeatures;
    std::vector<std::string> auauTightBDTCentINDControlFeatures;
    std::vector<std::string> auauTightBDTCentAsFeatFeatures;
    std::vector<std::string> auauTightBDTCentDepFeatures;
    std::string pre = "reference";
    std::string tight = "reference";
    std::string nonTight = "reference";
    std::string viewId;
    std::string materialize = "physics";
    std::string yamlText;
  };

  std::string defaultViewId(const Config& cfg);

  Config loadConfigFromPath(const std::string& path)
  {
    Config cfg;
    std::ifstream in(path);
    if (!in) return cfg;

    std::ostringstream text;
    text << in.rdbuf();
    cfg.yamlText = text.str();

    std::istringstream ss(cfg.yamlText);
    for (std::string line; std::getline(ss, line); )
    {
      line = trim(line);
      if (line.empty() || line[0] == '#') continue;
      if (startsWithKey(line, "use_vz_cut")) parseBool(afterColon(line), cfg.useVzCut);
      else if (startsWithKey(line, "vz_cut_cm")) parseDouble(afterColon(line), cfg.vzCutCm);
      else if (startsWithKey(line, "jes3_photon_pt_bins")) parseDoubles(afterColon(line), cfg.ptBins);
      else if (startsWithKey(line, "unfold_reco_photon_pt_bins")) parseDoubles(afterColon(line), cfg.unfoldRecoPhotonPtBins);
      else if (startsWithKey(line, "unfold_truth_photon_pt_bins")) parseDoubles(afterColon(line), cfg.unfoldTruthPhotonPtBins);
      else if (startsWithKey(line, "unfold_xj_bins")) parseDoubles(afterColon(line), cfg.unfoldXJBins);
      else if (startsWithKey(line, "jet_radii")) parseDoubles(afterColon(line), cfg.jetRadii);
      else if (startsWithKey(line, "centrality_edges")) parseInts(afterColon(line), cfg.centEdges);
      else if (startsWithKey(line, "jet_pt_min")) parseDouble(afterColon(line), cfg.minJetPt);
      else if (startsWithKey(line, "back_to_back_dphi_min_pi_fraction"))
      {
        double frac = 0.0;
        if (parseDouble(afterColon(line), frac)) cfg.minBackToBack = frac * M_PI;
      }
      else if (startsWithKey(line, "coneR")) parseDouble(afterColon(line), cfg.coneR);
      else if (startsWithKey(line, "isSlidingIso")) parseBool(afterColon(line), cfg.slidingIso);
      else if (startsWithKey(line, "fixedGeV")) parseDouble(afterColon(line), cfg.isoFixed);
      else if (startsWithKey(line, "isolation_wp"))
      {
        parseInlineMapDouble(line, "aGeV", cfg.isoA);
        parseInlineMapDouble(line, "bPerGeV", cfg.isoB);
        parseInlineMapDouble(line, "sideGapGeV", cfg.isoGap);
      }
      else if (startsWithKey(line, "A")) parseDouble(afterColon(line), cfg.isoA);
      else if (startsWithKey(line, "B")) parseDouble(afterColon(line), cfg.isoB);
      else if (startsWithKey(line, "sideband_gap")) parseDouble(afterColon(line), cfg.isoGap);
      else if (startsWithKey(line, "npb_cut")) parseDouble(afterColon(line), cfg.npbCut);
      else if (startsWithKey(line, "auau_npb_cut")) parseDouble(afterColon(line), cfg.auauNpbCut);
      else if (startsWithKey(line, "tight_bdt_min_intercept")) parseDouble(afterColon(line), cfg.tightBDTMinIntercept);
      else if (startsWithKey(line, "tight_bdt_min_slope")) parseDouble(afterColon(line), cfg.tightBDTMinSlope);
      else if (startsWithKey(line, "nontight_bdt_min_intercept")) parseDouble(afterColon(line), cfg.nonTightBDTMinIntercept);
      else if (startsWithKey(line, "nontight_bdt_min_slope")) parseDouble(afterColon(line), cfg.nonTightBDTMinSlope);
      else if (startsWithKey(line, "nontight_bdt_max_intercept")) parseDouble(afterColon(line), cfg.nonTightBDTMaxIntercept);
      else if (startsWithKey(line, "nontight_bdt_max_slope")) parseDouble(afterColon(line), cfg.nonTightBDTMaxSlope);
      else if (startsWithKey(line, "auau_tight_bdt_min_intercept")) parseDouble(afterColon(line), cfg.auauTightBDTMinIntercept);
      else if (startsWithKey(line, "auau_tight_bdt_min_slope")) parseDouble(afterColon(line), cfg.auauTightBDTMinSlope);
      else if (startsWithKey(line, "auau_tight_bdt_max")) parseDouble(afterColon(line), cfg.auauTightBDTMax);
      else if (startsWithKey(line, "auau_nontight_bdt_min_intercept")) parseDouble(afterColon(line), cfg.auauNonTightBDTMinIntercept);
      else if (startsWithKey(line, "auau_nontight_bdt_min_slope")) parseDouble(afterColon(line), cfg.auauNonTightBDTMinSlope);
      else if (startsWithKey(line, "auau_nontight_bdt_max_intercept")) parseDouble(afterColon(line), cfg.auauNonTightBDTMaxIntercept);
      else if (startsWithKey(line, "auau_nontight_bdt_max_slope")) parseDouble(afterColon(line), cfg.auauNonTightBDTMaxSlope);
      else if (startsWithKey(line, "auau_npb_model_file")) cfg.auauNpbModelFile = afterColon(line);
      else if (startsWithKey(line, "auau_npb_features")) parseStrings(afterColon(line), cfg.auauNpbFeatures);
      else if (startsWithKey(line, "auau_tight_bdt_model_file")) cfg.auauTightBDTModelFile = afterColon(line);
      else if (startsWithKey(line, "auau_tight_bdt_centINDcontrol_model_file")) cfg.auauTightBDTCentINDControlModelFile = afterColon(line);
      else if (startsWithKey(line, "auau_tight_bdt_centAsFeat_model_file")) cfg.auauTightBDTCentAsFeatModelFile = afterColon(line);
      else if (startsWithKey(line, "auau_tight_bdt_centDep_model_files")) parseStrings(afterColon(line), cfg.auauTightBDTCentDepModelFiles);
      else if (startsWithKey(line, "auau_tight_bdt_features")) parseStrings(afterColon(line), cfg.auauTightBDTFeatures);
      else if (startsWithKey(line, "auau_tight_bdt_centINDcontrol_features")) parseStrings(afterColon(line), cfg.auauTightBDTCentINDControlFeatures);
      else if (startsWithKey(line, "auau_tight_bdt_centAsFeat_features")) parseStrings(afterColon(line), cfg.auauTightBDTCentAsFeatFeatures);
      else if (startsWithKey(line, "auau_tight_bdt_centDep_features")) parseStrings(afterColon(line), cfg.auauTightBDTCentDepFeatures);
      else if (startsWithKey(line, "preselection")) cfg.pre = normalizePre(afterColon(line));
      else if (startsWithKey(line, "tight")) cfg.tight = normalizeTight(afterColon(line));
      else if (startsWithKey(line, "nonTight")) cfg.nonTight = normalizeNonTight(afterColon(line));
    }
    return cfg;
  }

  Config loadConfig()
  {
    const char* pathEnv = std::getenv("RJ_CONFIG_YAML");
    const std::string path = pathEnv && *pathEnv ? pathEnv : "macros/analysis_config.yaml";
    return loadConfigFromPath(path);
  }

  bool isInlineViewSpec(const std::string& s)
  {
    return s.rfind("INLINE_VIEW_V1", 0) == 0;
  }

  void applyInlineViewSpec(Config& cfg, const std::string& spec)
  {
    if (!isInlineViewSpec(spec)) return;
    std::stringstream ss(spec);
    for (std::string tok; std::getline(ss, tok, ';'); )
    {
      tok = trim(tok);
      if (tok.empty() || tok == "INLINE_VIEW_V1") continue;
      const std::size_t eq = tok.find('=');
      if (eq == std::string::npos) continue;
      const std::string key = trim(tok.substr(0, eq));
      const std::string val = trim(tok.substr(eq + 1));

      if (key == "jet_pt_min") parseDouble(val, cfg.minJetPt);
      else if (key == "back_to_back_dphi_min_pi_fraction")
      {
        double frac = 0.0;
        if (parseDouble(val, frac)) cfg.minBackToBack = frac * M_PI;
      }
      else if (key == "vz_cut_cm") parseDouble(val, cfg.vzCutCm);
      else if (key == "coneR") parseDouble(val, cfg.coneR);
      else if (key == "isSlidingIso") parseBool(val, cfg.slidingIso);
      else if (key == "fixedGeV") parseDouble(val, cfg.isoFixed);
      else if (key == "view_id") cfg.viewId = val;
    }
  }

  struct Entry
  {
    std::string outRoot;
    std::string cfgTag;
    std::string yamlPath;
    std::string pre;
    std::string tight;
    std::string nonTight;
    std::string viewId;
    std::string materialize = "physics";
    Config cfg;
  };

  std::vector<Entry> loadFanout(const std::string& outRoot, const Config& cfg)
  {
    std::vector<Entry> out;
    const char* pathEnv = std::getenv("RJ_REPLAY_FANOUT_FILE");
    if (!pathEnv || !*pathEnv) pathEnv = std::getenv("RJ_ID_FANOUT_FILE");
    if (!pathEnv || !*pathEnv)
    {
      Entry e;
      e.outRoot = outRoot;
      e.pre = cfg.pre;
      e.tight = cfg.tight;
      e.nonTight = cfg.nonTight;
      e.viewId = cfg.viewId.empty() ? defaultViewId(cfg) : cfg.viewId;
      e.materialize = cfg.materialize;
      e.cfg = cfg;
      out.push_back(e);
      return out;
    }

    std::ifstream in(pathEnv);
    for (std::string line; std::getline(in, line); )
    {
      line = trim(line);
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> cols;
      std::stringstream ss(line);
      for (std::string tok; std::getline(ss, tok, '|'); ) cols.push_back(trim(tok));
      if (cols.size() != 5 && cols.size() < 6) continue;

      Entry e;
      e.outRoot = cols[0];
      e.cfgTag = cols[1];
      e.cfg = cfg;
      if (cols.size() >= 6)
      {
        const std::string viewSpecOrYaml = cols[2];
        if (isInlineViewSpec(viewSpecOrYaml))
        {
          applyInlineViewSpec(e.cfg, viewSpecOrYaml);
        }
        else
        {
          e.yamlPath = viewSpecOrYaml;
          if (!e.yamlPath.empty()) e.cfg = loadConfigFromPath(e.yamlPath);
        }
        e.pre = normalizePre(cols[3]);
        e.tight = normalizeTight(cols[4]);
        e.nonTight = normalizeNonTight(cols[5]);
      }
      else
      {
        e.pre = normalizePre(cols[2]);
        e.tight = normalizeTight(cols[3]);
        e.nonTight = normalizeNonTight(cols[4]);
      }
      e.cfg.pre = e.pre;
      e.cfg.tight = e.tight;
      e.cfg.nonTight = e.nonTight;
      if (cols.size() >= 7) e.viewId = trim(cols[6]);
      if (cols.size() >= 8 && !trim(cols[7]).empty()) e.materialize = trim(cols[7]);
      if (e.viewId.empty()) e.viewId = defaultViewId(e.cfg);
      e.cfg.viewId = e.viewId;
      e.cfg.materialize = e.materialize;
      out.push_back(e);
    }
    if (out.empty())
    {
      Entry e;
      e.outRoot = outRoot;
      e.pre = cfg.pre;
      e.tight = cfg.tight;
      e.nonTight = cfg.nonTight;
      e.viewId = cfg.viewId.empty() ? defaultViewId(cfg) : cfg.viewId;
      e.materialize = cfg.materialize;
      e.cfg = cfg;
      out.push_back(e);
    }
    return out;
  }

  const std::vector<Entry>* gForcedFanoutEntries = nullptr;

  std::size_t replayMaxOpenOutputs()
  {
    const char* env = std::getenv("RJ_REPLAY_MAX_OPEN_OUTPUTS");
    if (!env || !*env) return 48;
    char* endptr = nullptr;
    const long v = std::strtol(env, &endptr, 10);
    return (endptr != env && v > 0) ? static_cast<std::size_t>(v) : 48;
  }

  bool envBool(const char* key, bool fallback = false)
  {
    const char* raw = std::getenv(key);
    if (!raw || !*raw) return fallback;
    std::string v = trim(raw);
    std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return std::tolower(c); });
    return v == "1" || v == "true" || v == "yes" || v == "on";
  }

  bool validateBranches(TChain& chain,
                        const std::string& treeName,
                        const std::vector<std::string>& required,
                        bool requiredTree = true)
  {
    if (chain.GetEntries() <= 0)
    {
      if (requiredTree)
      {
        std::cerr << "[RecoilJetsPoolReplay] Missing or empty pool tree: " << treeName << "\n";
        return false;
      }
      return true;
    }

    chain.LoadTree(0);
    bool ok = true;
    for (const std::string& br : required)
    {
      if (!chain.GetBranch(br.c_str()))
      {
        std::cerr << "[RecoilJetsPoolReplay] Pool schema mismatch: tree " << treeName
                  << " is missing required branch '" << br << "'\n";
        ok = false;
      }
    }
    return ok;
  }

  bool inOpen(double x, double lo, double hi)
  {
    return std::isfinite(x) && x > lo && x < hi;
  }

  int findPtBin(double pt, const std::vector<double>& bins)
  {
    for (std::size_t i = 0; i + 1 < bins.size(); ++i)
      if (pt >= bins[i] && pt < bins[i + 1]) return static_cast<int>(i);
    return -1;
  }

  int findCentBin(int cent, const std::vector<int>& edges)
  {
    for (std::size_t i = 0; i + 1 < edges.size(); ++i)
      if (cent >= edges[i] && cent < edges[i + 1]) return static_cast<int>(i);
    return -1;
  }

  long long scopedEventKey(const TChain& chain, long long rawEventKey)
  {
    std::uint64_t x = static_cast<std::uint64_t>(rawEventKey);
    const std::uint64_t y = static_cast<std::uint64_t>(chain.GetTreeNumber() + 1);
    x ^= y + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
    return static_cast<long long>(x);
  }

  std::string suffix(int ptIdx, int centIdx, bool isAuAu, const Config& cfg)
  {
    std::ostringstream os;
    if (ptIdx >= 0 && ptIdx + 1 < static_cast<int>(cfg.ptBins.size()))
      os << "_pT_" << std::fixed << std::setprecision(0) << cfg.ptBins[ptIdx] << "_" << cfg.ptBins[ptIdx + 1];
    if (isAuAu && centIdx >= 0 && centIdx + 1 < static_cast<int>(cfg.centEdges.size()))
      os << "_cent_" << cfg.centEdges[centIdx] << "_" << cfg.centEdges[centIdx + 1];
    return os.str();
  }

  bool passPre(const Entry& e, const Config& cfg, double pt, double weta, double et1, double e11e33, double e32e35, double npb, double auauNpb)
  {
    (void) pt;
    const bool ref = (e11e33 < PRE_E11E33_MAX) &&
                     inOpen(et1, PRE_ET1_MIN, PRE_ET1_MAX) &&
                     inOpen(e32e35, PRE_E32E35_MIN, PRE_E32E35_MAX) &&
                     (std::isfinite(weta) && weta < PRE_WETA_MAX);
    if (e.pre == "noPreCriteria") return true;
    if (e.pre == "onlyNPB") return std::isfinite(npb) && npb > cfg.npbCut;
    if (e.pre == "auauOnlyNPB") return std::isfinite(auauNpb) && auauNpb > cfg.auauNpbCut;
    if (e.pre == "newPPG12")
      return std::isfinite(npb) && npb > cfg.npbCut &&
             inOpen(et1, PRE_ET1_MIN, PRE_ET1_MAX) &&
             e11e33 < PRE_E11E33_MAX &&
             inOpen(e32e35, PRE_E32E35_MIN, PRE_E32E35_MAX);
    if (e.pre == "refPlusNPB") return ref && std::isfinite(npb) && npb > cfg.npbCut;
    return ref;
  }

  int tightTag(const Entry& e, const Config& cfg, double pt, double weta, double wphi, double et1, double e11e33, double e32e35, double tightBDT, double auauBDT)
  {
    if (e.tight == "newPPG12")
    {
      const double minT = cfg.tightBDTMinIntercept + cfg.tightBDTMinSlope * pt;
      if (std::isfinite(tightBDT) && tightBDT > minT) return 1;
      if (e.nonTight == "newPPG12")
      {
        const double lo = cfg.nonTightBDTMinIntercept + cfg.nonTightBDTMinSlope * pt;
        const double hi = cfg.nonTightBDTMaxIntercept + cfg.nonTightBDTMaxSlope * pt;
        return (std::isfinite(tightBDT) && tightBDT > lo && tightBDT < hi) ? 2 : 3;
      }
      return 2;
    }

    if (e.tight == "auauEmbeddedBDT" || e.tight == "centINDcontrol" || e.tight == "centAsFeat" || e.tight == "centDepBDTs")
    {
      const double loT = cfg.auauTightBDTMinIntercept + cfg.auauTightBDTMinSlope * pt;
      const bool isT = std::isfinite(auauBDT) && auauBDT > loT && auauBDT < cfg.auauTightBDTMax;
      if (isT) return 1;
      if (e.nonTight == "auauBDTSideband")
      {
        const double lo = cfg.auauNonTightBDTMinIntercept + cfg.auauNonTightBDTMinSlope * pt;
        const double hi = cfg.auauNonTightBDTMaxIntercept + cfg.auauNonTightBDTMaxSlope * pt;
        return (std::isfinite(auauBDT) && auauBDT > lo && auauBDT < hi) ? 2 : 3;
      }
      return 2;
    }

    const double wHi = 0.15 + 0.006 * pt;
    const int fails =
      (!inOpen(weta, TIGHT_W_LO, wHi)) +
      (!inOpen(wphi, TIGHT_W_LO, wHi)) +
      (!inOpen(e11e33, TIGHT_E11E33_MIN, TIGHT_E11E33_MAX)) +
      (!inOpen(et1, TIGHT_ET1_MIN, TIGHT_ET1_MAX)) +
      (!inOpen(e32e35, TIGHT_E32E35_MIN, TIGHT_E32E35_MAX));
    if (fails == 0) return 1;
    if (fails >= 2) return 2;
    return 3;
  }

  struct EventRec
  {
    int run = 0;
    int evt = 0;
    bool isSim = false;
    bool isAuAu = false;
    int centBin = -1;
    float centPercent = std::numeric_limits<float>::quiet_NaN();
    float vz = std::numeric_limits<float>::quiet_NaN();
    float weight = 1.0f;
    std::vector<std::string> triggers;
  };

  struct PhotonRec
  {
    long long eventKey = 0;
    int index = -1;
    float pt = std::numeric_limits<float>::quiet_NaN();
    float eta = std::numeric_limits<float>::quiet_NaN();
    float phi = std::numeric_limits<float>::quiet_NaN();
    float energy = std::numeric_limits<float>::quiet_NaN();
    float eiso = std::numeric_limits<float>::quiet_NaN();
    float eisoR03 = std::numeric_limits<float>::quiet_NaN();
    float eisoR04 = std::numeric_limits<float>::quiet_NaN();
    float iso03Emcal = std::numeric_limits<float>::quiet_NaN();
    float iso03Hcalin = std::numeric_limits<float>::quiet_NaN();
    float iso03Hcalout = std::numeric_limits<float>::quiet_NaN();
    float iso04Emcal = std::numeric_limits<float>::quiet_NaN();
    float iso04Hcalin = std::numeric_limits<float>::quiet_NaN();
    float iso04Hcalout = std::numeric_limits<float>::quiet_NaN();
    float weta = std::numeric_limits<float>::quiet_NaN();
    float wphi = std::numeric_limits<float>::quiet_NaN();
    float weta33 = std::numeric_limits<float>::quiet_NaN();
    float wphi33 = std::numeric_limits<float>::quiet_NaN();
    float weta35 = std::numeric_limits<float>::quiet_NaN();
    float wphi53 = std::numeric_limits<float>::quiet_NaN();
    float et1 = std::numeric_limits<float>::quiet_NaN();
    float et2 = std::numeric_limits<float>::quiet_NaN();
    float et3 = std::numeric_limits<float>::quiet_NaN();
    float et4 = std::numeric_limits<float>::quiet_NaN();
    float e11e33 = std::numeric_limits<float>::quiet_NaN();
    float e32e35 = std::numeric_limits<float>::quiet_NaN();
    float e11e22 = std::numeric_limits<float>::quiet_NaN();
    float e11e13 = std::numeric_limits<float>::quiet_NaN();
    float e11e15 = std::numeric_limits<float>::quiet_NaN();
    float e11e17 = std::numeric_limits<float>::quiet_NaN();
    float e11e31 = std::numeric_limits<float>::quiet_NaN();
    float e11e51 = std::numeric_limits<float>::quiet_NaN();
    float e11e71 = std::numeric_limits<float>::quiet_NaN();
    float e22e33 = std::numeric_limits<float>::quiet_NaN();
    float e22e35 = std::numeric_limits<float>::quiet_NaN();
    float e22e37 = std::numeric_limits<float>::quiet_NaN();
    float e22e53 = std::numeric_limits<float>::quiet_NaN();
    float w32 = std::numeric_limits<float>::quiet_NaN();
    float w52 = std::numeric_limits<float>::quiet_NaN();
    float w72 = std::numeric_limits<float>::quiet_NaN();
    float meanTime = std::numeric_limits<float>::quiet_NaN();
    float npb = std::numeric_limits<float>::quiet_NaN();
    float tightBDT = std::numeric_limits<float>::quiet_NaN();
    float auauNpb = std::numeric_limits<float>::quiet_NaN();
    float auauTightBDT = std::numeric_limits<float>::quiet_NaN();
    float mbdTime = std::numeric_limits<float>::quiet_NaN();
    float clusterMbdDeltaT = std::numeric_limits<float>::quiet_NaN();
    int npbHasAwayJet = 0;
    int npbLabel = -1;
    int isNPB = -1;
    int truthSignal = 0;
    int truthTrackId = -1;
    float truthPt = std::numeric_limits<float>::quiet_NaN();
    float truthEta = std::numeric_limits<float>::quiet_NaN();
    float truthPhi = std::numeric_limits<float>::quiet_NaN();
    std::vector<float> extraFeatures;
  };

  const std::vector<std::string>& analysisPoolExtraFeatureNames()
  {
    static const std::vector<std::string> names = {
      "cluster_energy", "cluster_et", "cluster_pt", "cluster_eta", "cluster_phi", "cluster_r", "cluster_z",
      "cluster_ecore", "cluster_prob", "cluster_chi2", "cluster_ntowers", "cluster_sum_tower_e",
      "cluster_max_tower_e", "cluster_max_tower_ieta", "cluster_max_tower_iphi",
      "e11", "e22", "e33", "e55", "e77", "e13", "e15", "e17", "e31", "e51", "e71",
      "e35", "e37", "e53", "e73", "e57", "e75", "e32", "e52", "e72",
      "et1", "et2", "et3", "et4",
      "e11_over_e33", "e32_over_e35", "e11_over_e22", "e11_over_e13", "e11_over_e15",
      "e11_over_e17", "e11_over_e31", "e11_over_e51", "e11_over_e71",
      "e22_over_e33", "e22_over_e35", "e22_over_e37", "e22_over_e53",
      "e11_over_e55", "e11_over_e77", "e33_over_e55", "e33_over_e77",
      "weta", "wphi", "weta_cog", "wphi_cog", "weta_cogx", "wphi_cogx",
      "weta33_cogx", "wphi33_cogx", "weta35_cogx", "wphi53_cogx",
      "w32", "w52", "w72", "detamax", "dphimax", "detacog", "dphicog", "drad",
      "ppg12_iso_axis_eta", "ppg12_iso_axis_phi", "vertex_z",
      "mean_time", "mbd_time", "cluster_mbd_delta_t",
      "npb_has_away_jet", "npb_label", "is_npb",
      "npb_score", "tight_bdt_score", "auau_npb_score", "auau_tight_bdt_score",
      "eiso", "eiso_r03", "eiso_r04",
      "iso_04_emcal", "iso_04_hcalin", "iso_04_hcalout",
      "iso_03_emcal", "iso_03_hcalin", "iso_03_hcalout",
      "iso_02_emcal", "iso_01_emcal", "iso_005_emcal",
      "ihcal_et", "ohcal_et", "ihcal_et22", "ohcal_et22", "ihcal_et33", "ohcal_et33",
      "ihcal_ieta", "ihcal_iphi", "ohcal_ieta", "ohcal_iphi",
      "cluster_tower_x_raw", "cluster_tower_y_raw", "cluster_tower_x_corr", "cluster_tower_y_corr",
      "cluster_incidence_alpha_phi", "cluster_incidence_alpha_eta"
    };
    return names;
  }

  std::string scoreCacheKey(const std::string& modelFile,
                            const std::vector<std::string>& features,
                            const EventRec& ev,
                            const PhotonRec& ph)
  {
    std::ostringstream os;
    os << modelFile << '|';
    for (const std::string& f : features) os << f << ',';
    os << '|' << ev.run << '|' << ev.evt << '|' << ph.index << '|'
       << std::setprecision(8) << ph.pt << '|' << ph.eta << '|' << ph.phi;
    return os.str();
  }

  double photonFeatureValue(const std::string& feature, const EventRec& ev, const PhotonRec& ph)
  {
    if (feature == "cluster_Et" || feature == "cluster_et" || feature == "ET") return ph.pt;
    if (feature == "cluster_Eta" || feature == "cluster_eta") return ph.eta;
    if (feature == "cluster_Phi" || feature == "cluster_phi") return ph.phi;
    if (feature == "vertexz" || feature == "vertex_z" || feature == "zvtx") return ev.vz;
    if (feature == "centrality" || feature == "cent") return std::isfinite(ev.centPercent) ? ev.centPercent : static_cast<float>(ev.centBin);
    if (feature == "event_weight") return ev.weight;
    if (feature == "reco_eiso") return std::isfinite(ph.eisoR04) ? ph.eisoR04 : ph.eiso;
    if (feature == "cluster_weta_cogx") return ph.weta;
    if (feature == "cluster_wphi_cogx") return ph.wphi;
    if (feature == "cluster_weta33_cogx") return ph.weta33;
    if (feature == "cluster_wphi33_cogx") return ph.wphi33;
    if (feature == "cluster_weta35_cogx") return ph.weta35;
    if (feature == "cluster_wphi53_cogx") return ph.wphi53;
    if (feature == "cluster_et1") return ph.et1;
    if (feature == "cluster_et2") return ph.et2;
    if (feature == "cluster_et3") return ph.et3;
    if (feature == "cluster_et4") return ph.et4;
    if (feature == "e11_over_e33") return ph.e11e33;
    if (feature == "e32_over_e35") return ph.e32e35;
    if (feature == "e11_over_e22") return ph.e11e22;
    if (feature == "e11_over_e13") return ph.e11e13;
    if (feature == "e11_over_e15") return ph.e11e15;
    if (feature == "e11_over_e17") return ph.e11e17;
    if (feature == "e11_over_e31") return ph.e11e31;
    if (feature == "e11_over_e51") return ph.e11e51;
    if (feature == "e11_over_e71") return ph.e11e71;
    if (feature == "e22_over_e33") return ph.e22e33;
    if (feature == "e22_over_e35") return ph.e22e35;
    if (feature == "e22_over_e37") return ph.e22e37;
    if (feature == "e22_over_e53") return ph.e22e53;
    if (feature == "cluster_w32") return ph.w32;
    if (feature == "cluster_w52") return ph.w52;
    if (feature == "cluster_w72") return ph.w72;
    if (feature == "mean_time" || feature == "cluster_mean_time") return ph.meanTime;
    if (feature == "mbd_time") return ph.mbdTime;
    if (feature == "cluster_mbd_delta_t") return ph.clusterMbdDeltaT;
    if (feature == "npb_score") return ph.npb;
    if (feature == "auau_npb_score") return ph.auauNpb;
    if (feature == "auau_tight_bdt_score") return ph.auauTightBDT;
    if (feature == "npb_has_away_jet") return ph.npbHasAwayJet;
    if (feature == "npb_label") return ph.npbLabel;
    if (feature == "is_npb") return ph.isNPB;
    const auto& extraNames = analysisPoolExtraFeatureNames();
    const auto it = std::find(extraNames.begin(), extraNames.end(), feature);
    if (it != extraNames.end())
    {
      const std::size_t idx = static_cast<std::size_t>(std::distance(extraNames.begin(), it));
      if (idx < ph.extraFeatures.size()) return ph.extraFeatures[idx];
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  struct ReplayBDTModels
  {
    std::unordered_map<std::string, std::unique_ptr<TMVA::Experimental::RBDT>> models;
    std::unordered_map<std::string, double> scores;
    std::unordered_map<std::string, int> warned;

    TMVA::Experimental::RBDT* model(const std::string& file)
    {
      if (file.empty()) return nullptr;
      auto it = models.find(file);
      if (it != models.end()) return it->second.get();
      try
      {
        auto ptr = std::make_unique<TMVA::Experimental::RBDT>("myBDT", file);
        auto* raw = ptr.get();
        models.emplace(file, std::move(ptr));
        std::cout << "[RecoilJetsPoolReplay] loaded replay BDT model " << file << "\n";
        return raw;
      }
      catch (const std::exception& e)
      {
        if (!warned[file]++)
          std::cerr << "[RecoilJetsPoolReplay] failed to load replay BDT model " << file
                    << ": " << e.what() << "\n";
      }
      catch (...)
      {
        if (!warned[file]++)
          std::cerr << "[RecoilJetsPoolReplay] failed to load replay BDT model " << file << "\n";
      }
      return nullptr;
    }

    double eval(const std::string& file,
                const std::vector<std::string>& features,
                const EventRec& ev,
                const PhotonRec& ph)
    {
      if (file.empty() || features.empty()) return std::numeric_limits<double>::quiet_NaN();
      const std::string key = scoreCacheKey(file, features, ev, ph);
      if (auto it = scores.find(key); it != scores.end()) return it->second;
      auto* m = model(file);
      if (!m) return std::numeric_limits<double>::quiet_NaN();
      std::vector<float> x;
      x.reserve(features.size());
      for (const std::string& f : features)
      {
        const double v = photonFeatureValue(f, ev, ph);
        if (!std::isfinite(v))
        {
          std::ostringstream msg;
          msg << "[RecoilJetsPoolReplay] missing feature for replay BDT"
              << " feature=" << f
              << " model=" << file
              << " run=" << ev.run
              << " evt=" << ev.evt
              << " photon_index=" << ph.index
              << " pt=" << ph.pt;
          throw std::runtime_error(msg.str());
        }
        x.push_back(static_cast<float>(v));
      }
      const auto y = m->Compute(x);
      scores[key] = (!y.empty() && std::isfinite(y[0])) ? static_cast<double>(y[0])
                                                        : std::numeric_limits<double>::quiet_NaN();
      return scores[key];
    }
  };

  std::vector<std::string> fallbackFeatures(const std::vector<std::string>& preferred,
                                            const std::vector<std::string>& fallback)
  {
    return preferred.empty() ? fallback : preferred;
  }

  double replayAuAuNpbScore(ReplayBDTModels& bdt, const Entry& e, const EventRec& ev, const PhotonRec& ph)
  {
    const Config& cfg = e.cfg;
    if (!cfg.auauNpbModelFile.empty() && !cfg.auauNpbFeatures.empty())
    {
      const double s = bdt.eval(cfg.auauNpbModelFile, cfg.auauNpbFeatures, ev, ph);
      if (std::isfinite(s)) return s;
    }
    return ph.auauNpb;
  }

  double replayAuAuTightScore(ReplayBDTModels& bdt, const Entry& e, const EventRec& ev, const PhotonRec& ph)
  {
    const Config& cfg = e.cfg;
    std::string modelFile;
    std::vector<std::string> features;

    if (e.tight == "centDepBDTs")
    {
      features = fallbackFeatures(cfg.auauTightBDTCentDepFeatures, cfg.auauTightBDTFeatures);
      if (cfg.auauTightBDTCentDepModelFiles.size() + 1 == cfg.centEdges.size())
      {
        const double cent = std::isfinite(ev.centPercent) ? ev.centPercent : static_cast<double>(ev.centBin);
        for (std::size_t i = 0; i + 1 < cfg.centEdges.size(); ++i)
        {
          const bool inBin = cent >= cfg.centEdges[i] &&
            (cent < cfg.centEdges[i + 1] || (i + 2 == cfg.centEdges.size() && cent <= cfg.centEdges[i + 1]));
          if (inBin)
          {
            modelFile = cfg.auauTightBDTCentDepModelFiles[i];
            break;
          }
        }
      }
    }
    else if (e.tight == "centINDcontrol")
    {
      modelFile = cfg.auauTightBDTCentINDControlModelFile.empty() ? cfg.auauTightBDTModelFile
                                                                  : cfg.auauTightBDTCentINDControlModelFile;
      features = fallbackFeatures(cfg.auauTightBDTCentINDControlFeatures, cfg.auauTightBDTFeatures);
    }
    else if (e.tight == "centAsFeat")
    {
      modelFile = cfg.auauTightBDTCentAsFeatModelFile.empty() ? cfg.auauTightBDTModelFile
                                                              : cfg.auauTightBDTCentAsFeatModelFile;
      features = withCentralityFeature(fallbackFeatures(cfg.auauTightBDTCentAsFeatFeatures, cfg.auauTightBDTFeatures));
    }
    else if (e.tight == "auauEmbeddedBDT")
    {
      modelFile = cfg.auauTightBDTModelFile;
      features = cfg.auauTightBDTFeatures;
    }

    if (!modelFile.empty() && !features.empty())
    {
      const double s = bdt.eval(modelFile, features, ev, ph);
      if (std::isfinite(s)) return s;
    }
    return ph.auauTightBDT;
  }

  struct JetRec
  {
    long long eventKey = 0;
    std::string rKey;
    int isTruth = 0;
    float pt = std::numeric_limits<float>::quiet_NaN();
    float rawPt = std::numeric_limits<float>::quiet_NaN();
    float areaSubPt = std::numeric_limits<float>::quiet_NaN();
    float eta = std::numeric_limits<float>::quiet_NaN();
    float phi = std::numeric_limits<float>::quiet_NaN();
    float mass = std::numeric_limits<float>::quiet_NaN();
    float area = std::numeric_limits<float>::quiet_NaN();
  };

  struct TruthPhotonRec
  {
    long long eventKey = 0;
    int trackId = -1;
    int barcode = -1;
    float pt = std::numeric_limits<float>::quiet_NaN();
    float eta = std::numeric_limits<float>::quiet_NaN();
    float phi = std::numeric_limits<float>::quiet_NaN();
    float iso = std::numeric_limits<float>::quiet_NaN();
  };

  struct View
  {
    Entry entry;
    Config cfg;
    std::string viewId;
    std::string materialize = "physics";
    bool legacyTopLevel = false;
  };

  struct Output
  {
    std::string outRoot;
    TFile* file = nullptr;
    std::map<std::string, std::map<std::string, TH1*> > h;
    std::vector<View> views;
    const View* activeView = nullptr;
  };

  std::string centSuffix(int centIdx, bool isAuAu, const Config& cfg)
  {
    if (!isAuAu || centIdx < 0 || centIdx + 1 >= static_cast<int>(cfg.centEdges.size())) return "";
    std::ostringstream os;
    os << "_cent_" << cfg.centEdges[centIdx] << "_" << cfg.centEdges[centIdx + 1];
    return os.str();
  }

  std::string compactNumber(double v, int scale = 100)
  {
    if (scale > 0)
    {
      const double x = v * static_cast<double>(scale);
      if (std::fabs(x - std::round(x)) < 1e-6)
      {
        std::ostringstream os;
        os << static_cast<long long>(std::llround(x));
        return os.str();
      }
    }
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << v;
    std::string s = os.str();
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
    std::replace(s.begin(), s.end(), '.', 'p');
    std::replace(s.begin(), s.end(), '-', 'm');
    return s;
  }

  std::string defaultViewId(const Config& cfg)
  {
    std::ostringstream os;
    os << "jetPt" << compactNumber(cfg.minJetPt, 1)
       << "_dphi";
    const double frac = cfg.minBackToBack / M_PI;
    if (std::fabs(frac - 0.5) < 1e-6) os << "pi2";
    else if (std::fabs(frac - 0.875) < 1e-6) os << "7pi8";
    else os << "piFrac" << compactNumber(frac, 1000);
    os << "_vz" << compactNumber(cfg.vzCutCm, 1)
       << "_";
    if (cfg.slidingIso) os << "isoSliding";
    else os << "isoFixed" << compactNumber(cfg.isoFixed, 1);
    os << "_coneR" << compactNumber(cfg.coneR, 100);
    return os.str();
  }

  double photonIsoForView(const PhotonRec& ph, const Config& cfg)
  {
    if (std::fabs(cfg.coneR - 0.3) < 1e-6 && std::isfinite(ph.eisoR03)) return ph.eisoR03;
    if (std::fabs(cfg.coneR - 0.4) < 1e-6 && std::isfinite(ph.eisoR04)) return ph.eisoR04;
    return ph.eiso;
  }

  double photonIsoPartForView(const PhotonRec& ph, const Config& cfg, const std::string& part)
  {
    if (std::fabs(cfg.coneR - 0.3) < 1e-6)
    {
      if (part == "emcal") return ph.iso03Emcal;
      if (part == "hcalin") return ph.iso03Hcalin;
      if (part == "hcalout") return ph.iso03Hcalout;
    }
    if (std::fabs(cfg.coneR - 0.4) < 1e-6)
    {
      if (part == "emcal") return ph.iso04Emcal;
      if (part == "hcalin") return ph.iso04Hcalin;
      if (part == "hcalout") return ph.iso04Hcalout;
    }
    return std::numeric_limits<double>::quiet_NaN();
  }

  double jetRFromKey(const std::string& rKey)
  {
    if (rKey.size() >= 2 && (rKey[0] == 'r' || rKey[0] == 'R'))
    {
      char* endptr = nullptr;
      const std::string tail = rKey.substr(1);
      const double value = std::strtod(tail.c_str(), &endptr);
      if (endptr != tail.c_str()) return value / 100.0;
    }
    return 0.4;
  }

  double jetEtaAbsMaxForRKey(const std::string& rKey)
  {
    return std::max(0.0, 1.1 - jetRFromKey(rKey));
  }

  double dR(double eta1, double phi1, double eta2, double phi2)
  {
    const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
    const double deta = eta1 - eta2;
    return std::sqrt(deta * deta + dphi * dphi);
  }

  TDirectory* ensureDir(Output& out, const std::string& trig)
  {
    if (!out.file)
    {
      return gDirectory;
    }
    out.file->cd();
    TDirectory* dir = out.file->GetDirectory(trig.c_str());
    if (!dir) dir = out.file->mkdir(trig.c_str());
    if (!dir) return nullptr;
    const View* view = out.activeView;
    if (!view || view->legacyTopLevel || view->viewId.empty())
    {
      dir->cd();
      return dir;
    }
    TDirectory* viewsDir = dir->GetDirectory("views");
    if (!viewsDir) viewsDir = dir->mkdir("views");
    if (!viewsDir)
    {
      dir->cd();
      return dir;
    }
    TDirectory* viewDir = viewsDir->GetDirectory(view->viewId.c_str());
    if (!viewDir) viewDir = viewsDir->mkdir(view->viewId.c_str());
    if (viewDir) viewDir->cd();
    return viewDir ? viewDir : dir;
  }

  TDirectory* ensurePath(TFile& file, const std::string& path)
  {
    file.cd();
    TDirectory* dir = &file;
    std::stringstream ss(path);
    for (std::string part; std::getline(ss, part, '/'); )
    {
      part = trim(part);
      if (part.empty()) continue;
      TDirectory* next = dir->GetDirectory(part.c_str());
      if (!next) next = dir->mkdir(part.c_str());
      if (!next) return dir;
      dir = next;
    }
    dir->cd();
    return dir;
  }

  std::string histScopeKey(const Output& out, const std::string& trig)
  {
    const View* view = out.activeView;
    if (!view || view->legacyTopLevel || view->viewId.empty()) return trig;
    return trig + "/views/" + view->viewId;
  }

  TH1I* countHist(Output& out, const std::string& trig, const std::string& name)
  {
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1I*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1I(name.c_str(), (name + ";count;entries").c_str(), 1, 0.5, 1.5);
    h->GetXaxis()->SetBinLabel(1, "count");
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH1F* vertexHist(Output& out, const std::string& trig)
  {
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find("h_vertexZ"); it != m.end()) return dynamic_cast<TH1F*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1F("h_vertexZ", "h_vertexZ;v_{z} [cm];Entries", 240, -60.0, 60.0);
    h->Sumw2();
    m["h_vertexZ"] = h;
    return h;
  }

  TH1F* hist1(Output& out, const std::string& trig, const std::string& name,
              const std::string& title, int nb, double lo, double hi)
  {
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1F*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1F(name.c_str(), title.c_str(), nb, lo, hi);
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH1I* hist1I(Output& out, const std::string& trig, const std::string& name,
               const std::string& title, int nb, double lo, double hi)
  {
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1I*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1I(name.c_str(), title.c_str(), nb, lo, hi);
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH1F* hist1Var(Output& out, const std::string& trig, const std::string& name,
                 const std::string& title, const std::vector<double>& bins)
  {
    if (bins.size() < 2) return nullptr;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1F*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1F(name.c_str(), title.c_str(), static_cast<int>(bins.size()) - 1, bins.data());
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH2F* hist2(Output& out, const std::string& trig, const std::string& name,
              const std::string& title,
              const std::vector<double>& xbins, int ny, double ylo, double yhi)
  {
    if (xbins.size() < 2) return nullptr;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH2F*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH2F(name.c_str(), title.c_str(),
                       static_cast<int>(xbins.size()) - 1, xbins.data(),
                       ny, ylo, yhi);
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH2F* hist2Var(Output& out, const std::string& trig, const std::string& name,
                 const std::string& title,
                 const std::vector<double>& xbins, const std::vector<double>& ybins)
  {
    if (xbins.size() < 2 || ybins.size() < 2) return nullptr;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH2F*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH2F(name.c_str(), title.c_str(),
                       static_cast<int>(xbins.size()) - 1, xbins.data(),
                       static_cast<int>(ybins.size()) - 1, ybins.data());
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH2F* hist2Uniform(Output& out, const std::string& trig, const std::string& name,
                     const std::string& title,
                     int nx, double xlo, double xhi,
                     int ny, double ylo, double yhi)
  {
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH2F*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH2F(name.c_str(), title.c_str(), nx, xlo, xhi, ny, ylo, yhi);
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH3F* hist3VarUniform(Output& out, const std::string& trig, const std::string& name,
                        const std::string& title,
                        const std::vector<double>& xbins,
                        int ny, double ylo, double yhi,
                        int nz, double zlo, double zhi)
  {
    if (xbins.size() < 2) return nullptr;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH3F*>(it->second);
    std::vector<double> ybins(ny + 1), zbins(nz + 1);
    for (int i = 0; i <= ny; ++i) ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);
    for (int i = 0; i <= nz; ++i) zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);
    ensureDir(out, trig);
    auto* h = new TH3F(name.c_str(), title.c_str(),
                       static_cast<int>(xbins.size()) - 1, xbins.data(),
                       ny, ybins.data(),
                       nz, zbins.data());
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TProfile* profileVar(Output& out, const std::string& trig, const std::string& name,
                       const std::string& title, const std::vector<double>& xbins,
                       double ylo, double yhi)
  {
    if (xbins.size() < 2) return nullptr;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TProfile*>(it->second);
    ensureDir(out, trig);
    auto* h = new TProfile(name.c_str(), title.c_str(),
                           static_cast<int>(xbins.size()) - 1, xbins.data(),
                           ylo, yhi);
    m[name] = h;
    return h;
  }

  TH1I* abcdHist(Output& out, const std::string& trig, const std::string& base, const std::string& sfx)
  {
    const std::string name = base + sfx;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1I*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1I(name.c_str(), (name + ";count;entries").c_str(), 1, 0.5, 1.5);
    h->GetXaxis()->SetBinLabel(1, "count");
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH1F* isoHist(Output& out, const std::string& trig, const std::string& base, const std::string& sfx,
                const std::string& xAxis = "E_{T}^{iso} [GeV]")
  {
    return hist1(out, trig, base + sfx, base + sfx + ";" + xAxis + ";Entries", 170, -5.0, 12.0);
  }

  TH1I* isoDecisionHist(Output& out, const std::string& trig, const std::string& sfx)
  {
    const std::string name = "h_isoDecision" + sfx;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1I*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1I(name.c_str(), (name + ";Isolation decision;Entries").c_str(), 2, 0.5, 2.5);
    h->GetXaxis()->SetBinLabel(1, "PASS");
    h->GetXaxis()->SetBinLabel(2, "FAIL");
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH1F* ptGammaABCDHist(Output& out, const std::string& trig, const std::string& base, int centIdx,
                        bool isAuAu, const Config& cfg)
  {
    const std::string name = base + suffix(-1, centIdx, isAuAu, cfg);
    return hist1Var(out, trig, name, name + ";p_{T}^{#gamma} [GeV];Entries", cfg.ptBins);
  }

  TH1F* ssHist(Output& out, const std::string& trig, const std::string& varKey,
               const std::string& tagKey, const std::string& sfx)
  {
    int nb = 120;
    double lo = 0.0, hi = 1.2;
    if (varKey == "npbScore" || varKey == "tightBDTScore" || varKey == "auauNpbScore" || varKey == "auauTightBDTScore")
    {
      lo = -1.0;
      hi = 1.2;
    }
    else if (varKey == "meanTime" || varKey == "clusterMbdDeltaT")
    {
      lo = -50.0;
      hi = 50.0;
      nb = 200;
    }
    const std::string name = "h_ss_" + varKey + "_" + tagKey + sfx;
    return hist1(out, trig, name, name + ";" + varKey + ";Entries", nb, lo, hi);
  }

  TH3F* hist3Legacy(Output& out, const std::string& trig, const std::string& name,
                    const std::string& title,
                    const std::vector<double>& xbins,
                    int ny, double ylo, double yhi,
                    int nz, double zlo, double zhi)
  {
    return hist3VarUniform(out, trig, name, title, xbins, ny, ylo, yhi, nz, zlo, zhi);
  }

  TProfile3D* profile3Legacy(Output& out, const std::string& trig, const std::string& name,
                             const std::string& title,
                             const std::vector<double>& xbins,
                             int ny, double ylo, double yhi,
                             int nz, double zlo, double zhi)
  {
    if (xbins.size() < 2) return nullptr;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TProfile3D*>(it->second);
    std::vector<double> ybins(ny + 1), zbins(nz + 1);
    for (int i = 0; i <= ny; ++i) ybins[i] = ylo + (yhi - ylo) * (static_cast<double>(i) / ny);
    for (int i = 0; i <= nz; ++i) zbins[i] = zlo + (zhi - zlo) * (static_cast<double>(i) / nz);
    ensureDir(out, trig);
    auto* p = new TProfile3D(name.c_str(), title.c_str(),
                             static_cast<int>(xbins.size()) - 1, xbins.data(),
                             ny, ybins.data(),
                             nz, zbins.data(),
                             "");
    m[name] = p;
    return p;
  }

  TH1I* sigABCDHist(Output& out, const std::string& trig, const std::string& sfx)
  {
    const std::string name = "h_sigABCD_MC" + sfx;
    auto& m = out.h[histScopeKey(out, trig)];
    if (auto it = m.find(name); it != m.end()) return dynamic_cast<TH1I*>(it->second);
    ensureDir(out, trig);
    auto* h = new TH1I(name.c_str(), (name + ";Reco ABCD region for matched truth-signal photons;Entries").c_str(), 4, 0.5, 4.5);
    h->GetXaxis()->SetBinLabel(1, "A");
    h->GetXaxis()->SetBinLabel(2, "B");
    h->GetXaxis()->SetBinLabel(3, "C");
    h->GetXaxis()->SetBinLabel(4, "D");
    h->Sumw2();
    m[name] = h;
    return h;
  }

  TH1F* xJHist(Output& out, const std::string& trig, const std::string& rKey, const std::string& sfx)
  {
    return hist1(out, trig, "h_xJ_" + rKey + sfx,
                 "h_xJ_" + rKey + sfx + ";x_{J};Entries", 60, 0.0, 3.0);
  }

  TH1F* jet1Hist(Output& out, const std::string& trig, const std::string& rKey, const std::string& sfx)
  {
    return hist1(out, trig, "h_jet1Pt_" + rKey + sfx,
                 "h_jet1Pt_" + rKey + sfx + ";p_{T}^{jet1} [GeV];Entries", 120, 0.0, 60.0);
  }

  TH1F* jet2Hist(Output& out, const std::string& trig, const std::string& rKey, const std::string& sfx)
  {
    return hist1(out, trig, "h_jet2Pt_" + rKey + sfx,
                 "h_jet2Pt_" + rKey + sfx + ";p_{T}^{jet2} [GeV];Entries", 120, 0.0, 60.0);
  }

  TH1F* alphaHist(Output& out, const std::string& trig, const std::string& rKey, const std::string& sfx)
  {
    return hist1(out, trig, "h_alpha_" + rKey + sfx,
                 "h_alpha_" + rKey + sfx + ";#alpha=p_{T}^{jet2}/p_{T}^{#gamma};Entries", 60, 0.0, 2.0);
  }

  int globalBin2D(const std::vector<double>& xbins, const std::vector<double>& ybins, double x, double y)
  {
    if (xbins.size() < 2 || ybins.size() < 2) return -1;
    TH2F tmp("rjpool_tmp_globalbin", "", static_cast<int>(xbins.size()) - 1, xbins.data(),
             static_cast<int>(ybins.size()) - 1, ybins.data());
    return tmp.FindBin(x, y);
  }

  void fillPhotonCounts(Output& out, const EventRec& ev, int centIdx, const Config& cfg,
                        const PhotonRec* leadRecoTight, const TruthPhotonRec* leadTruth)
  {
    for (const std::string& trig : ev.triggers)
    {
      const std::string cs = centSuffix(centIdx, ev.isAuAu, cfg);
      const bool haveRecoForNominal = leadRecoTight && (!ev.isSim || leadRecoTight->truthSignal);
      const bool recoFake = ev.isSim && leadRecoTight && !leadRecoTight->truthSignal;
      const bool haveTruthMatchedReco = ev.isSim && leadRecoTight && leadRecoTight->truthSignal;
      if (haveRecoForNominal)
      {
        if (auto* h = hist1Var(out, trig, "h_unfoldRecoPho_pTgamma" + cs,
                               "h_unfoldRecoPho_pTgamma" + cs + ";p_{T}^{#gamma,reco} [GeV];N_{#gamma}^{reco}",
                               cfg.unfoldRecoPhotonPtBins))
          h->Fill(leadRecoTight->pt, ev.weight);
      }
      if (leadTruth)
      {
        if (auto* h = hist1Var(out, trig, "h_unfoldTruthPho_pTgamma" + cs,
                               "h_unfoldTruthPho_pTgamma" + cs + ";p_{T}^{#gamma,truth} [GeV];N_{#gamma}^{truth}",
                               cfg.unfoldTruthPhotonPtBins))
          h->Fill(leadTruth->pt, ev.weight);
      }
      if (leadTruth && haveTruthMatchedReco)
      {
        if (auto* h = hist2Var(out, trig, "h2_unfoldResponsePho_pTgamma" + cs,
                               "h2_unfoldResponsePho_pTgamma" + cs + ";p_{T}^{#gamma,truth} [GeV];p_{T}^{#gamma,reco} [GeV]",
                               cfg.unfoldTruthPhotonPtBins, cfg.unfoldRecoPhotonPtBins))
          h->Fill(leadTruth->pt, leadRecoTight->pt, ev.weight);
      }
      else if (recoFake)
      {
        if (auto* h = hist1Var(out, trig, "h_unfoldRecoPhoFakes_pTgamma" + cs,
                               "h_unfoldRecoPhoFakes_pTgamma" + cs + ";p_{T}^{#gamma,reco} [GeV];Reco fakes",
                               cfg.unfoldRecoPhotonPtBins))
          h->Fill(leadRecoTight->pt, ev.weight);
      }
      if (leadTruth && ev.isSim && !haveTruthMatchedReco)
      {
        if (auto* h = hist1Var(out, trig, "h_unfoldTruthPhoMisses_pTgamma" + cs,
                               "h_unfoldTruthPhoMisses_pTgamma" + cs + ";p_{T}^{#gamma,truth} [GeV];Truth misses",
                               cfg.unfoldTruthPhotonPtBins))
          h->Fill(leadTruth->pt, ev.weight);
      }
    }
  }

  double ssValue(const PhotonRec& ph, const std::string& key)
  {
    if (key == "weta") return ph.weta;
    if (key == "wphi") return ph.wphi;
    if (key == "et1") return ph.et1;
    if (key == "e11e33") return ph.e11e33;
    if (key == "e32e35") return ph.e32e35;
    if (key == "npbScore") return ph.npb;
    if (key == "tightBDTScore") return ph.tightBDT;
    if (key == "auauNpbScore") return ph.auauNpb;
    if (key == "auauTightBDTScore") return ph.auauTightBDT;
    if (key == "meanTime") return ph.meanTime;
    if (key == "clusterMbdDeltaT") return ph.clusterMbdDeltaT;
    return std::numeric_limits<double>::quiet_NaN();
  }

  void fillSSSet(Output& out, const std::string& trig, const PhotonRec& ph,
                 const std::string& tag, const std::string& sfx, double weight)
  {
    static const std::vector<std::string> vars = {
      "weta", "wphi", "et1", "e11e33", "e32e35",
      "npbScore", "tightBDTScore", "auauNpbScore", "auauTightBDTScore",
      "meanTime", "clusterMbdDeltaT"
    };
    for (const std::string& v : vars)
    {
      const double x = ssValue(ph, v);
      if (!std::isfinite(x)) continue;
      if (auto* h = ssHist(out, trig, v, tag, sfx)) h->Fill(x, weight);
    }
  }

  void fillPhotonLegacyQA(Output& out, const EventRec& ev, int centIdx, const Config& cfg,
                          const PhotonRec& ph, int tightTagValue,
                          bool iso, bool nonIso, int regionBin)
  {
    const int ptIdx = findPtBin(ph.pt, cfg.ptBins);
    if (ptIdx < 0) return;

    const std::string sfx = suffix(ptIdx, centIdx, ev.isAuAu, cfg);
    const double isoValue = photonIsoForView(ph, cfg);
    const double isoEmcal = photonIsoPartForView(ph, cfg, "emcal");
    const double isoHcalin = photonIsoPartForView(ph, cfg, "hcalin");
    const double isoHcalout = photonIsoPartForView(ph, cfg, "hcalout");
    const bool validIso = std::isfinite(isoValue) && isoValue < 1e8;

    for (const std::string& trig : ev.triggers)
    {
      if (validIso)
      {
        if (auto* h = isoHist(out, trig, "h_Eiso", sfx)) h->Fill(isoValue, ev.weight);
        if (auto* h = isoDecisionHist(out, trig, sfx)) h->Fill(iso ? 1 : 2, ev.weight);
      }
      if (std::isfinite(isoEmcal))
        if (auto* h = isoHist(out, trig, "h_Eiso_emcal", sfx, "E_{T}^{iso,EMCal} [GeV]")) h->Fill(isoEmcal, ev.weight);
      if (std::isfinite(isoHcalin))
        if (auto* h = isoHist(out, trig, "h_Eiso_hcalin", sfx, "E_{T}^{iso,IHCAL} [GeV]")) h->Fill(isoHcalin, ev.weight);
      if (std::isfinite(isoHcalout))
        if (auto* h = isoHist(out, trig, "h_Eiso_hcalout", sfx, "E_{T}^{iso,OHCAL} [GeV]")) h->Fill(isoHcalout, ev.weight);

      fillSSSet(out, trig, ph, "inclusive", sfx, ev.weight);
      if (iso) fillSSSet(out, trig, ph, "iso", sfx, ev.weight);
      else if (nonIso) fillSSSet(out, trig, ph, "nonIso", sfx, ev.weight);

      const std::string mcSuffix = ev.isSim ? (ph.truthSignal ? "_sig" : "_bkg") : "";
      fillSSSet(out, trig, ph, "pre" + mcSuffix, sfx, ev.weight);

      if (tightTagValue == 1 || tightTagValue == 2)
      {
        const bool tight = (tightTagValue == 1);
        const std::string tag = tight ? "tight" : "nonTight";
        const std::string eisoTag = tight ? "h_Eiso_tight" : "h_Eiso_nonTight";
        const std::string eisoEmTag = tight ? "h_Eiso_emcal_tight" : "h_Eiso_emcal_nonTight";
        const std::string eisoHiTag = tight ? "h_Eiso_hcalin_tight" : "h_Eiso_hcalin_nonTight";
        const std::string eisoHoTag = tight ? "h_Eiso_hcalout_tight" : "h_Eiso_hcalout_nonTight";
        if (validIso)
          if (auto* h = isoHist(out, trig, eisoTag, sfx)) h->Fill(isoValue, ev.weight);
        if (std::isfinite(isoEmcal))
          if (auto* h = isoHist(out, trig, eisoEmTag, sfx, "E_{T}^{iso,EMCal} [GeV]")) h->Fill(isoEmcal, ev.weight);
        if (std::isfinite(isoHcalin))
          if (auto* h = isoHist(out, trig, eisoHiTag, sfx, "E_{T}^{iso,IHCAL} [GeV]")) h->Fill(isoHcalin, ev.weight);
        if (std::isfinite(isoHcalout))
          if (auto* h = isoHist(out, trig, eisoHoTag, sfx, "E_{T}^{iso,OHCAL} [GeV]")) h->Fill(isoHcalout, ev.weight);
        fillSSSet(out, trig, ph, tag + mcSuffix, sfx, ev.weight);
      }

      if (ph.truthSignal && validIso)
      {
        if (auto* h = isoHist(out, trig, "h_EisoReco_truthSigMatched", sfx)) h->Fill(isoValue, ev.weight);
        if (tightTagValue == 1)
          if (auto* h = isoHist(out, trig, "h_EisoReco_truthSigMatched_tight", sfx)) h->Fill(isoValue, ev.weight);
      }

      if (regionBin > 0)
      {
        const char region = " ABCD"[regionBin];
        const std::string regionTag =
          (regionBin == 1) ? "isIsolated_isTight" :
          (regionBin == 2) ? "notIsolated_isTight" :
          (regionBin == 3) ? "isIsolated_notTight" : "notIsolated_notTight";
        fillSSSet(out, trig, ph, regionTag, sfx, ev.weight);

        const std::string hIsoBase = std::string("h_Eiso_ABCD_") + region;
        if (validIso)
          if (auto* h = isoHist(out, trig, hIsoBase, sfx)) h->Fill(isoValue, ev.weight);

        const std::string hPtBase = std::string("h_pTgamma_ABCD_") + region;
        if (auto* h = ptGammaABCDHist(out, trig, hPtBase, centIdx, ev.isAuAu, cfg)) h->Fill(ph.pt, ev.weight);
      }
    }
  }

  void fillInclusiveJetQA(Output& out, const EventRec& ev, int centIdx, const Config& cfg,
                          const std::vector<JetRec>& jets)
  {
    const std::string cs = centSuffix(centIdx, ev.isAuAu, cfg);
    std::map<std::string, std::vector<const JetRec*>> byR;
    for (const JetRec& j : jets)
    {
      if (j.isTruth || j.rKey.empty()) continue;
      byR[j.rKey].push_back(&j);
    }

    for (const auto& [rKey, rJets] : byR)
    {
      const double etaAbsMax = jetEtaAbsMaxForRKey(rKey);
      int nJetsFid = 0;
      double ht = 0.0;
      const JetRec* lead = nullptr;
      const JetRec* sub = nullptr;

      for (const JetRec* jp : rJets)
      {
        if (!jp) continue;
        const JetRec& j = *jp;
        if (!std::isfinite(j.pt) || !std::isfinite(j.eta) || !std::isfinite(j.phi)) continue;
        if (j.pt < cfg.minJetPt) continue;

        for (const std::string& trig : ev.triggers)
        {
          if (auto* h = hist1(out, trig, "h_jetPt_all_" + rKey + cs,
                              "h_jetPt_all_" + rKey + cs + ";p_{T}^{jet} [GeV];Entries", 200, 0.0, 100.0))
            h->Fill(j.pt, ev.weight);
          if (auto* h = hist1(out, trig, "h_jetEta_all_" + rKey + cs,
                              "h_jetEta_all_" + rKey + cs + ";#eta^{jet};Entries", 120, -1.2, 1.2))
            h->Fill(j.eta, ev.weight);
          if (auto* h = hist1(out, trig, "h_jetPhi_all_" + rKey + cs,
                              "h_jetPhi_all_" + rKey + cs + ";#phi^{jet} [rad];Entries", 128, -M_PI, M_PI))
            h->Fill(TVector2::Phi_mpi_pi(j.phi), ev.weight);
          if (auto* h = hist2Uniform(out, trig, "h_jetEtaPhi_all_" + rKey + cs,
                                     "h_jetEtaPhi_all_" + rKey + cs + ";#eta^{jet};#phi^{jet} [rad]",
                                     60, -1.2, 1.2, 64, -M_PI, M_PI))
            h->Fill(j.eta, TVector2::Phi_mpi_pi(j.phi), ev.weight);
          if (std::isfinite(j.mass) && j.mass >= 0.0)
          {
            if (auto* h = hist1(out, trig, "h_jetMass_all_" + rKey + cs,
                                "h_jetMass_all_" + rKey + cs + ";m^{jet} [GeV];Entries", 120, 0.0, 30.0))
              h->Fill(j.mass, ev.weight);
            if (auto* h = hist2Uniform(out, trig, "h_jetMassVsPt_all_" + rKey + cs,
                                       "h_jetMassVsPt_all_" + rKey + cs + ";p_{T}^{jet} [GeV];m^{jet} [GeV]",
                                       100, 0.0, 100.0, 80, 0.0, 20.0))
              h->Fill(j.pt, j.mass, ev.weight);
          }
        }

        if (std::fabs(j.eta) >= etaAbsMax) continue;
        ++nJetsFid;
        ht += j.pt;
        if (!lead || j.pt > lead->pt)
        {
          sub = lead;
          lead = &j;
        }
        else if (!sub || j.pt > sub->pt)
        {
          sub = &j;
        }

        for (const std::string& trig : ev.triggers)
        {
          if (auto* h = hist1(out, trig, "h_jetPt_incl_" + rKey + cs,
                              "h_jetPt_incl_" + rKey + cs + ";p_{T}^{jet} [GeV];Entries", 200, 0.0, 100.0))
            h->Fill(j.pt, ev.weight);
          if (auto* h = hist1(out, trig, "h_jetEta_incl_" + rKey + cs,
                              "h_jetEta_incl_" + rKey + cs + ";#eta^{jet};Entries", 120, -1.2, 1.2))
            h->Fill(j.eta, ev.weight);
          if (auto* h = hist1(out, trig, "h_jetPhi_incl_" + rKey + cs,
                              "h_jetPhi_incl_" + rKey + cs + ";#phi^{jet} [rad];Entries", 128, -M_PI, M_PI))
            h->Fill(TVector2::Phi_mpi_pi(j.phi), ev.weight);
          if (auto* h = hist2Uniform(out, trig, "h_jetEtaPhi_incl_" + rKey + cs,
                                     "h_jetEtaPhi_incl_" + rKey + cs + ";#eta^{jet};#phi^{jet} [rad]",
                                     60, -1.2, 1.2, 64, -M_PI, M_PI))
            h->Fill(j.eta, TVector2::Phi_mpi_pi(j.phi), ev.weight);
          if (std::isfinite(j.mass) && j.mass >= 0.0)
          {
            if (auto* h = hist1(out, trig, "h_jetMass_incl_" + rKey + cs,
                                "h_jetMass_incl_" + rKey + cs + ";m^{jet} [GeV];Entries", 120, 0.0, 30.0))
              h->Fill(j.mass, ev.weight);
            if (auto* h = hist2Uniform(out, trig, "h_jetMassVsPt_incl_" + rKey + cs,
                                       "h_jetMassVsPt_incl_" + rKey + cs + ";p_{T}^{jet} [GeV];m^{jet} [GeV]",
                                       100, 0.0, 100.0, 80, 0.0, 20.0))
              h->Fill(j.pt, j.mass, ev.weight);
          }
        }
      }

      for (const std::string& trig : ev.triggers)
      {
        if (auto* h = hist1I(out, trig, "h_nJets_" + rKey + cs,
                             "h_nJets_" + rKey + cs + ";N_{jets};Entries", 100, 0.0, 100.0))
          h->Fill(nJetsFid, ev.weight);
        if (auto* h = hist1(out, trig, "h_HT_" + rKey + cs,
                            "h_HT_" + rKey + cs + ";H_{T} [GeV];Entries", 300, 0.0, 150.0))
          h->Fill(ht, ev.weight);
        if (lead)
        {
          if (auto* h = hist1(out, trig, "h_leadJetPt_" + rKey + cs,
                              "h_leadJetPt_" + rKey + cs + ";p_{T}^{lead jet} [GeV];Entries", 200, 0.0, 100.0))
            h->Fill(lead->pt, ev.weight);
          if (auto* h = hist1(out, trig, "h_leadJetEta_" + rKey + cs,
                              "h_leadJetEta_" + rKey + cs + ";#eta^{lead jet};Entries", 120, -1.2, 1.2))
            h->Fill(lead->eta, ev.weight);
          if (auto* h = hist1(out, trig, "h_leadJetPhi_" + rKey + cs,
                              "h_leadJetPhi_" + rKey + cs + ";#phi^{lead jet} [rad];Entries", 128, -M_PI, M_PI))
            h->Fill(TVector2::Phi_mpi_pi(lead->phi), ev.weight);
        }
        if (sub)
        {
          if (auto* h = hist1(out, trig, "h_subleadJetPt_" + rKey + cs,
                              "h_subleadJetPt_" + rKey + cs + ";p_{T}^{sublead jet} [GeV];Entries", 200, 0.0, 100.0))
            h->Fill(sub->pt, ev.weight);
        }
      }
    }
  }

  void fillRecoJetReplay(Output& out, const EventRec& ev, int centIdx, const Config& cfg,
                         const PhotonRec& leadPho, const std::vector<JetRec>& jets, bool sidebandC)
  {
    const int ptIdx = findPtBin(leadPho.pt, cfg.ptBins);
    if (ptIdx < 0) return;
    const std::string ptCentSfx = suffix(ptIdx, centIdx, ev.isAuAu, cfg);
    const std::string cs = centSuffix(centIdx, ev.isAuAu, cfg);

    std::map<std::string, std::vector<JetRec>> byR;
    for (const JetRec& j : jets)
      if (!j.isTruth) byR[j.rKey].push_back(j);

    for (const auto& [rKey, rJets] : byR)
    {
      const double etaAbsMax = jetEtaAbsMaxForRKey(rKey);
      int nPassPt = 0, nPassEta = 0, nPassDphi = 0;
      double all1Pt = -1.0;
      double all2Pt = -1.0;
      const JetRec* all1 = nullptr;
      double maxDphi = -1.0;

      std::vector<const JetRec*> fidJets;
      std::vector<char> fidIsRecoil;

      for (const JetRec& j : rJets)
      {
        if (!std::isfinite(j.pt) || !std::isfinite(j.eta) || !std::isfinite(j.phi)) continue;
        if (j.pt < cfg.minJetPt) continue;
        ++nPassPt;

        if (dR(j.eta, j.phi, leadPho.eta, leadPho.phi) < 0.4) continue;

        if (j.pt > all1Pt)
        {
          all2Pt = all1Pt;
          all1Pt = j.pt;
          all1 = &j;
        }
        else if (j.pt > all2Pt)
        {
          all2Pt = j.pt;
        }

        if (std::fabs(j.eta) < etaAbsMax)
        {
          ++nPassEta;
          const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(j.phi - leadPho.phi));
          if (std::isfinite(dphiAbs) && dphiAbs > maxDphi) maxDphi = dphiAbs;
          const bool isRecoil = dphiAbs >= cfg.minBackToBack;
          if (isRecoil) ++nPassDphi;
          fidJets.push_back(&j);
          fidIsRecoil.push_back(isRecoil ? 1 : 0);
        }
      }

      const JetRec* recoil1 = nullptr;
      if (all1 && std::fabs(all1->eta) < etaAbsMax)
      {
        const double dphiLead = std::fabs(TVector2::Phi_mpi_pi(all1->phi - leadPho.phi));
        if (std::isfinite(dphiLead) && dphiLead >= cfg.minBackToBack) recoil1 = all1;
      }

      if (!sidebandC)
      {
        int status = 0;
        if (nPassPt == 0) status = 1;
        else if (nPassEta == 0) status = 2;
        else if (nPassDphi == 0 || !recoil1) status = 3;
        else status = 4;

        for (const std::string& trig : ev.triggers)
        {
          if (auto* h = hist2(out, trig, "h_match_status_vs_pTgamma_" + rKey + cs,
                              "h_match_status_vs_pTgamma_" + rKey + cs + ";p_{T}^{#gamma} [GeV];match status",
                              cfg.ptBins, 4, 0.5, 4.5))
            h->Fill(leadPho.pt, status, ev.weight);
          if (auto* h = hist2(out, trig, "h_match_maxdphi_vs_pTgamma_" + rKey + cs,
                              "h_match_maxdphi_vs_pTgamma_" + rKey + cs + ";p_{T}^{#gamma} [GeV];max |#Delta#phi(#gamma, jet)| [rad]",
                              cfg.ptBins, 64, 0.0, M_PI))
            h->Fill(leadPho.pt, maxDphi >= 0.0 ? maxDphi : -0.01, ev.weight);
          if (auto* p = profileVar(out, trig, "p_nRecoilJets_vs_pTgamma_" + rKey + cs,
                                   "p_nRecoilJets_vs_pTgamma_" + rKey + cs + ";p_{T}^{#gamma} [GeV];#LT N_{recoil jets} #GT",
                                   cfg.ptBins, 0.0, 50.0))
            p->Fill(leadPho.pt, static_cast<double>(nPassDphi), ev.weight);
        }
      }

      for (std::size_t i = 0; i < fidJets.size(); ++i)
      {
        if (!fidIsRecoil[i]) continue;
        const JetRec* j = fidJets[i];
        if (!j || !std::isfinite(j->pt) || j->pt <= 0.0) continue;
        const double xJ = j->pt / leadPho.pt;
        const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(j->phi - leadPho.phi));
        for (const std::string& trig : ev.triggers)
        {
          if (sidebandC)
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldReco_pTgamma_xJ_incl_sidebandC_" + rKey + cs,
                                   "h2_unfoldReco_pTgamma_xJ_incl_sidebandC_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco}",
                                   cfg.unfoldRecoPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadPho.pt, xJ, ev.weight);
          }
          else
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldReco_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldReco_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco}",
                                   cfg.unfoldRecoPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadPho.pt, xJ, ev.weight);
            if (auto* h = hist2(out, trig, "h2_unfoldReco_pTgamma_dphi_incl_" + rKey + cs,
                                "h2_unfoldReco_pTgamma_dphi_incl_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];|#Delta#phi(#gamma,jet)| [rad]",
                                cfg.unfoldRecoPhotonPtBins, 64, 0.0, M_PI))
              h->Fill(leadPho.pt, dphiAbs, ev.weight);
          }
        }
      }

      if (!sidebandC && recoil1)
      {
        double jet2Pt = 0.0;
        const JetRec* jet2 = nullptr;
        const double overlapR = std::max(0.30, jetRFromKey(rKey));
        for (const JetRec* j : fidJets)
        {
          if (!j || j == recoil1) continue;
          if (dR(j->eta, j->phi, leadPho.eta, leadPho.phi) <= overlapR) continue;
          if (j->pt > jet2Pt)
          {
            jet2Pt = j->pt;
            jet2 = j;
          }
        }
        const double xJ = recoil1->pt / leadPho.pt;
        const double alpha = jet2Pt > 0.0 ? jet2Pt / leadPho.pt : 0.0;
        const double dphiSel = std::fabs(TVector2::Phi_mpi_pi(recoil1->phi - leadPho.phi));

        for (const std::string& trig : ev.triggers)
        {
          if (auto* h = hist2(out, trig, "h_match_dphi_vs_pTgamma_" + rKey + cs,
                              "h_match_dphi_vs_pTgamma_" + rKey + cs + ";p_{T}^{#gamma} [GeV];|#Delta#phi(#gamma,recoil jet)| [rad]",
                              cfg.ptBins, 64, 0.0, M_PI))
            h->Fill(leadPho.pt, dphiSel, ev.weight);
          if (auto* h = hist2(out, trig, "h_recoilIsLeading_vs_pTgamma_" + rKey + cs,
                              "h_recoilIsLeading_vs_pTgamma_" + rKey + cs + ";p_{T}^{#gamma} [GeV];recoil jet is leading jet",
                              cfg.ptBins, 2, -0.5, 1.5))
            h->Fill(leadPho.pt, (all1 && recoil1 == all1) ? 1.0 : 0.0, ev.weight);
          if (auto* h = xJHist(out, trig, rKey, ptCentSfx)) h->Fill(xJ, ev.weight);
          if (auto* h = jet1Hist(out, trig, rKey, ptCentSfx)) h->Fill(recoil1->pt, ev.weight);
          if (auto* h = jet2Hist(out, trig, rKey, ptCentSfx)) h->Fill(jet2Pt, ev.weight);
          if (auto* h = alphaHist(out, trig, rKey, ptCentSfx)) h->Fill(alpha, ev.weight);
          if (auto* h = hist1(out, trig, "h_jet1Eta_sel_" + rKey + ptCentSfx,
                              "h_jet1Eta_sel_" + rKey + ptCentSfx + ";#eta^{jet1};Entries", 120, -1.2, 1.2))
            h->Fill(recoil1->eta, ev.weight);
          if (auto* h = hist1(out, trig, "h_jet1Phi_sel_" + rKey + ptCentSfx,
                              "h_jet1Phi_sel_" + rKey + ptCentSfx + ";#phi^{jet1} [rad];Entries", 128, -M_PI, M_PI))
            h->Fill(TVector2::Phi_mpi_pi(recoil1->phi), ev.weight);
          if (std::isfinite(recoil1->mass) && recoil1->mass >= 0.0)
            if (auto* h = hist1(out, trig, "h_jet1Mass_sel_" + rKey + ptCentSfx,
                                "h_jet1Mass_sel_" + rKey + ptCentSfx + ";m^{jet1} [GeV];Entries", 120, 0.0, 30.0))
              h->Fill(recoil1->mass, ev.weight);
          if (jet2)
          {
            if (auto* h = hist1(out, trig, "h_jet2Eta_sel_" + rKey + ptCentSfx,
                                "h_jet2Eta_sel_" + rKey + ptCentSfx + ";#eta^{jet2};Entries", 120, -1.2, 1.2))
              h->Fill(jet2->eta, ev.weight);
            if (auto* h = hist1(out, trig, "h_jet2Phi_sel_" + rKey + ptCentSfx,
                                "h_jet2Phi_sel_" + rKey + ptCentSfx + ";#phi^{jet2} [rad];Entries", 128, -M_PI, M_PI))
              h->Fill(TVector2::Phi_mpi_pi(jet2->phi), ev.weight);
          }
          if (auto* h = hist3VarUniform(out, trig, "h_JES3_pT_xJ_alpha_" + rKey + cs,
                                        "h_JES3_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma} [GeV];x_{J};#alpha",
                                        cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
            h->Fill(leadPho.pt, xJ, alpha, ev.weight);
          if (auto* h = hist3VarUniform(out, trig, "h_JES3_pT_jet1Pt_alpha_" + rKey + cs,
                                        "h_JES3_pT_jet1Pt_alpha_" + rKey + cs + ";p_{T}^{#gamma} [GeV];p_{T}^{jet1} [GeV];#alpha",
                                        cfg.ptBins, 120, 0.0, 60.0, 40, 0.0, 2.0))
            h->Fill(leadPho.pt, recoil1->pt, alpha, ev.weight);
          if (auto* h = hist3Legacy(out, trig, "h_Jet13_pTgamma_eta_phi_recoilJet1_" + rKey,
                                    "h_Jet13_pTgamma_eta_phi_recoilJet1_" + rKey + ";p_{T}^{#gamma} [GeV];#eta^{jet1};#phi^{jet1} [rad]",
                                    cfg.ptBins, 60, -1.1, 1.1, 64, -M_PI, M_PI))
            h->Fill(leadPho.pt, recoil1->eta, TVector2::Phi_mpi_pi(recoil1->phi), ev.weight);
          if (auto* p = profile3Legacy(out, trig, "p_Balance3_pTgamma_eta_phi_" + rKey,
                                       "p_Balance3_pTgamma_eta_phi_" + rKey + ";p_{T}^{#gamma} [GeV];#eta^{jet1};#phi^{jet1} [rad];#LT x_{J} #GT",
                                       cfg.ptBins, 60, -1.1, 1.1, 64, -M_PI, M_PI))
            p->Fill(leadPho.pt, recoil1->eta, TVector2::Phi_mpi_pi(recoil1->phi), xJ, ev.weight);
          if (ev.isSim && leadPho.truthSignal)
          {
            if (auto* h = hist3VarUniform(out, trig, "h_JES3RecoTruthPhoTagged_pT_xJ_alpha_" + rKey + cs,
                                          "h_JES3RecoTruthPhoTagged_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J}^{reco};#alpha^{reco}",
                                          cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
              h->Fill(leadPho.pt, xJ, alpha, ev.weight);
          }
        }
      }
    }
  }

  void fillTruthJetReplay(Output& out, const EventRec& ev, int centIdx, const Config& cfg,
                          const TruthPhotonRec& truthPho, const std::vector<JetRec>& jets)
  {
    const std::string cs = centSuffix(centIdx, ev.isAuAu, cfg);
    std::map<std::string, std::vector<JetRec>> byR;
    for (const JetRec& j : jets)
      if (j.isTruth) byR[j.rKey].push_back(j);

    for (const auto& [rKey, truthJets] : byR)
    {
      const double etaAbsMax = jetEtaAbsMaxForRKey(rKey);
      const JetRec* leadTruth = nullptr;
      double leadTruthPt = -1.0;
      double subTruthPt = -1.0;

      for (const JetRec& j : truthJets)
      {
        if (!std::isfinite(j.pt) || !std::isfinite(j.eta) || !std::isfinite(j.phi)) continue;
        if (j.pt < cfg.minJetPt) continue;
        if (dR(j.eta, j.phi, truthPho.eta, truthPho.phi) < 0.4) continue;

        const bool fid = std::fabs(j.eta) < etaAbsMax;
        const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(j.phi - truthPho.phi));
        const bool recoil = fid && dphiAbs >= cfg.minBackToBack;

        if (recoil)
        {
          const double xJ = j.pt / truthPho.pt;
          for (const std::string& trig : ev.triggers)
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldTruth_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldTruth_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}",
                                   cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(truthPho.pt, xJ, ev.weight);
            if (auto* h = hist2(out, trig, "h2_unfoldTruth_pTgamma_dphi_incl_" + rKey + cs,
                                "h2_unfoldTruth_pTgamma_dphi_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];|#Delta#phi(#gamma^{truth},jet^{truth})| [rad]",
                                cfg.unfoldTruthPhotonPtBins, 64, 0.0, M_PI))
              h->Fill(truthPho.pt, dphiAbs, ev.weight);
          }
        }

        if (j.pt > leadTruthPt)
        {
          subTruthPt = leadTruthPt;
          leadTruthPt = j.pt;
          leadTruth = &j;
        }
        else if (j.pt > subTruthPt)
        {
          subTruthPt = j.pt;
        }
      }

      if (!leadTruth) continue;
      if (!(std::fabs(leadTruth->eta) < etaAbsMax)) continue;
      const double dphiLead = std::fabs(TVector2::Phi_mpi_pi(leadTruth->phi - truthPho.phi));
      if (!std::isfinite(dphiLead) || dphiLead < cfg.minBackToBack) continue;

      const double xJ = leadTruth->pt / truthPho.pt;
      const double alpha = subTruthPt > 0.0 ? subTruthPt / truthPho.pt : 0.0;
      for (const std::string& trig : ev.triggers)
      {
        if (auto* h = hist3VarUniform(out, trig, "h_JES3TruthPure_pT_xJ_alpha_" + rKey + cs,
                                      "h_JES3TruthPure_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J}^{truth};#alpha^{truth}",
                                      cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
          h->Fill(truthPho.pt, xJ, alpha, ev.weight);
        if (auto* h = hist3VarUniform(out, trig, "h_JES3TruthPure_pT_jet1Pt_alpha_" + rKey + cs,
                                      "h_JES3TruthPure_pT_jet1Pt_alpha_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];p_{T}^{jet1,truth} [GeV];#alpha^{truth}",
                                      cfg.ptBins, 120, 0.0, 60.0, 40, 0.0, 2.0))
          h->Fill(truthPho.pt, leadTruth->pt, alpha, ev.weight);
      }
    }
  }

  struct MatchedJetPair
  {
    const JetRec* reco = nullptr;
    const JetRec* truth = nullptr;
    double dr = std::numeric_limits<double>::quiet_NaN();
  };

  std::vector<MatchedJetPair> greedyMatchJets(const std::vector<const JetRec*>& reco,
                                              const std::vector<const JetRec*>& truth,
                                              double maxDR)
  {
    struct Cand { double dr; int ir; int it; };
    std::vector<Cand> cands;
    for (int ir = 0; ir < static_cast<int>(reco.size()); ++ir)
    {
      for (int it = 0; it < static_cast<int>(truth.size()); ++it)
      {
        const double dr = dR(reco[ir]->eta, reco[ir]->phi, truth[it]->eta, truth[it]->phi);
        if (std::isfinite(dr) && dr < maxDR) cands.push_back({dr, ir, it});
      }
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.dr < b.dr; });
    std::vector<char> usedReco(reco.size(), 0), usedTruth(truth.size(), 0);
    std::vector<MatchedJetPair> out;
    for (const Cand& c : cands)
    {
      if (usedReco[c.ir] || usedTruth[c.it]) continue;
      usedReco[c.ir] = 1;
      usedTruth[c.it] = 1;
      out.push_back({reco[c.ir], truth[c.it], c.dr});
    }
    return out;
  }

  void fillSimUnfoldResponseReplay(Output& out, const EventRec& ev, int centIdx, const Config& cfg,
                                   const PhotonRec* leadReco, const TruthPhotonRec* leadTruth,
                                   const std::vector<JetRec>& jets)
  {
    if (!ev.isSim || !leadTruth) return;
    const std::string cs = centSuffix(centIdx, ev.isAuAu, cfg);

    std::map<std::string, std::vector<const JetRec*>> recoByR;
    std::map<std::string, std::vector<const JetRec*>> truthByR;
    for (const JetRec& j : jets)
    {
      if (j.rKey.empty() || !std::isfinite(j.pt) || !std::isfinite(j.eta) || !std::isfinite(j.phi)) continue;
      if (j.isTruth) truthByR[j.rKey].push_back(&j);
      else recoByR[j.rKey].push_back(&j);
    }

    for (const auto& kv : truthByR)
    {
      const std::string& rKey = kv.first;
      const double etaAbsMax = jetEtaAbsMaxForRKey(rKey);
      std::vector<const JetRec*> truthRecoil;
      std::vector<const JetRec*> recoRecoil;

      for (const JetRec* tj : kv.second)
      {
        if (!tj || tj->pt < cfg.minJetPt) continue;
        if (std::fabs(tj->eta) >= etaAbsMax) continue;
        if (dR(tj->eta, tj->phi, leadTruth->eta, leadTruth->phi) < 0.4) continue;
        const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(tj->phi - leadTruth->phi));
        if (dphiAbs < cfg.minBackToBack) continue;
        truthRecoil.push_back(tj);
      }

      if (leadReco)
      {
        for (const JetRec* rj : recoByR[rKey])
        {
          if (!rj || rj->pt < cfg.minJetPt) continue;
          if (std::fabs(rj->eta) >= etaAbsMax) continue;
          if (dR(rj->eta, rj->phi, leadReco->eta, leadReco->phi) < 0.4) continue;
          const double dphiAbs = std::fabs(TVector2::Phi_mpi_pi(rj->phi - leadReco->phi));
          if (dphiAbs < cfg.minBackToBack) continue;
          recoRecoil.push_back(rj);
        }
      }

      auto pairs = greedyMatchJets(recoRecoil, truthRecoil, 0.3);
      std::vector<char> truthMatched(truthRecoil.size(), 0), recoMatched(recoRecoil.size(), 0);
      for (const MatchedJetPair& p : pairs)
      {
        auto ti = std::find(truthRecoil.begin(), truthRecoil.end(), p.truth);
        auto ri = std::find(recoRecoil.begin(), recoRecoil.end(), p.reco);
        if (ti != truthRecoil.end()) truthMatched[std::distance(truthRecoil.begin(), ti)] = 1;
        if (ri != recoRecoil.end()) recoMatched[std::distance(recoRecoil.begin(), ri)] = 1;
      }

      for (const std::string& trig : ev.triggers)
      {
        for (std::size_t it = 0; it < truthRecoil.size(); ++it)
        {
          const JetRec* tj = truthRecoil[it];
          const double xJt = tj->pt / leadTruth->pt;
          if (auto* h = hist2Var(out, trig, "h2_unfoldJetEffDen_pTgamma_xJ_incl_" + rKey + cs,
                                 "h2_unfoldJetEffDen_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}",
                                 cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins))
            h->Fill(leadTruth->pt, xJt, ev.weight);
          if (truthMatched[it])
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldJetEffNum_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldJetEffNum_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}",
                                   cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadTruth->pt, xJt, ev.weight);
            if (auto* h = hist2Var(out, trig, "h2_unfoldTruthMatched_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldTruthMatched_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}",
                                   cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadTruth->pt, xJt, ev.weight);
          }
          else
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldTruthMisses_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldTruthMisses_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}",
                                   cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadTruth->pt, xJt, ev.weight);
            if (auto* h = hist2Var(out, trig, "h2_unfoldTruthMisses_typeB_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldTruthMisses_typeB_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J#gamma}^{truth}",
                                   cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadTruth->pt, xJt, ev.weight);
          }
        }

        for (std::size_t ir = 0; ir < recoRecoil.size(); ++ir)
        {
          const JetRec* rj = recoRecoil[ir];
          const double xJr = leadReco ? rj->pt / leadReco->pt : std::numeric_limits<double>::quiet_NaN();
          if (!std::isfinite(xJr)) continue;
          if (recoMatched[ir])
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldRecoMatched_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldRecoMatched_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco}",
                                   cfg.unfoldRecoPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadReco->pt, xJr, ev.weight);
          }
          else
          {
            if (auto* h = hist2Var(out, trig, "h2_unfoldRecoFakes_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldRecoFakes_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco}",
                                   cfg.unfoldRecoPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadReco->pt, xJr, ev.weight);
            if (auto* h = hist2Var(out, trig, "h2_unfoldRecoFakes_typeB_pTgamma_xJ_incl_" + rKey + cs,
                                   "h2_unfoldRecoFakes_typeB_pTgamma_xJ_incl_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J#gamma}^{reco}",
                                   cfg.unfoldRecoPhotonPtBins, cfg.unfoldXJBins))
              h->Fill(leadReco->pt, xJr, ev.weight);
          }
        }

        for (const MatchedJetPair& p : pairs)
        {
          if (!p.reco || !p.truth || !leadReco) continue;
          const double xJr = p.reco->pt / leadReco->pt;
          const double xJt = p.truth->pt / leadTruth->pt;
          const int gTruth = globalBin2D(cfg.unfoldTruthPhotonPtBins, cfg.unfoldXJBins, leadTruth->pt, xJt);
          const int gReco = globalBin2D(cfg.unfoldRecoPhotonPtBins, cfg.unfoldXJBins, leadReco->pt, xJr);
          if (gTruth >= 0 && gReco >= 0)
          {
            if (auto* h = hist2Uniform(out, trig, "h2_unfoldResponse_pTgamma_xJ_incl_" + rKey + cs,
                                       "h2_unfoldResponse_pTgamma_xJ_incl_" + rKey + cs + ";global bin truth;global bin reco",
                                       (static_cast<int>(cfg.unfoldTruthPhotonPtBins.size()) + 1) * (static_cast<int>(cfg.unfoldXJBins.size()) + 1), -0.5,
                                       (static_cast<int>(cfg.unfoldTruthPhotonPtBins.size()) + 1) * (static_cast<int>(cfg.unfoldXJBins.size()) + 1) - 0.5,
                                       (static_cast<int>(cfg.unfoldRecoPhotonPtBins.size()) + 1) * (static_cast<int>(cfg.unfoldXJBins.size()) + 1), -0.5,
                                       (static_cast<int>(cfg.unfoldRecoPhotonPtBins.size()) + 1) * (static_cast<int>(cfg.unfoldXJBins.size()) + 1) - 0.5))
              h->Fill(gTruth, gReco, ev.weight);
          }
          if (auto* h = hist1(out, trig, "h_unfoldJetMatch_dR_" + rKey + cs,
                              "h_unfoldJetMatch_dR_" + rKey + cs + ";#DeltaR(reco jet, truth jet);Counts",
                              60, 0.0, 0.3))
            h->Fill(p.dr, ev.weight);
          if (p.truth->pt > 0.0)
          {
            if (auto* h = hist2Uniform(out, trig, "h2_unfoldJetPtResponse_pTtruth_ratio_" + rKey + cs,
                                       "h2_unfoldJetPtResponse_pTtruth_ratio_" + rKey + cs + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco}/p_{T}^{jet,truth}",
                                       120, 0.0, 60.0, 120, 0.0, 2.0))
              h->Fill(p.truth->pt, p.reco->pt / p.truth->pt, ev.weight);
            if (auto* h = hist2Uniform(out, trig, "h2_unfoldJetPtResponseAll_pTtruth_ratio_" + rKey + cs,
                                       "h2_unfoldJetPtResponseAll_pTtruth_ratio_" + rKey + cs + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco}/p_{T}^{jet,truth}",
                                       120, 0.0, 60.0, 120, 0.0, 2.0))
              h->Fill(p.truth->pt, p.reco->pt / p.truth->pt, ev.weight);
          }
        }

        const JetRec* truthLead = nullptr;
        for (const JetRec* tj : truthRecoil)
          if (!truthLead || tj->pt > truthLead->pt) truthLead = tj;
        const JetRec* recoLead = nullptr;
        for (const JetRec* rj : recoRecoil)
          if (!recoLead || rj->pt > recoLead->pt) recoLead = rj;
        if (truthLead)
        {
          if (auto* h = hist1Var(out, trig, "h_leadTruthRecoilMatch_den_pTgammaTruth_" + rKey + cs,
                                 "h_leadTruthRecoilMatch_den_pTgammaTruth_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];Den",
                                 cfg.ptBins))
            h->Fill(leadTruth->pt, ev.weight);
          const JetRec* truthLeadMatchedReco = nullptr;
          double truthLeadDR = std::numeric_limits<double>::infinity();
          for (const MatchedJetPair& p : pairs)
          {
            if (p.truth == truthLead)
            {
              truthLeadMatchedReco = p.reco;
              truthLeadDR = p.dr;
              break;
            }
          }
          const bool num = recoLead && truthLeadMatchedReco && recoLead == truthLeadMatchedReco;
          const bool missA = !num && truthLeadMatchedReco;
          const bool missB = !num && !truthLeadMatchedReco;
          const std::string cls = num ? "num" : (missA ? "missA" : "missB");
          if (auto* h = hist1Var(out, trig, "h_leadTruthRecoilMatch_" + cls + "_pTgammaTruth_" + rKey + cs,
                                 "h_leadTruthRecoilMatch_" + cls + "_pTgammaTruth_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];Entries",
                                 cfg.ptBins))
            h->Fill(leadTruth->pt, ev.weight);
          if (missA)
          {
            if (auto* h = hist1Var(out, trig, "h_leadTruthRecoilMatch_missA1_pTgammaTruth_" + rKey + cs,
                                   "h_leadTruthRecoilMatch_missA1_pTgammaTruth_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];Entries",
                                   cfg.ptBins))
              h->Fill(leadTruth->pt, ev.weight);
          }
          if (recoLead)
          {
            const double dphiReco = std::fabs(TVector2::Phi_mpi_pi(recoLead->phi - leadTruth->phi));
            const double xJReco = recoLead->pt / leadTruth->pt;
            if (auto* h = hist2(out, trig, "h2_leadTruthRecoilMatch_dphiRecoJet1_" + cls + "_pTgammaTruth_" + rKey + cs,
                                "h2_leadTruthRecoilMatch_dphiRecoJet1_" + cls + "_pTgammaTruth_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];|#Delta#phi|",
                                cfg.ptBins, 64, 0.0, M_PI))
              h->Fill(leadTruth->pt, dphiReco, ev.weight);
            if (auto* h = hist2(out, trig, "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_" + cls + "_pTgammaTruth_" + rKey + cs,
                                "h2_leadTruthRecoilMatch_dRRecoJet1_vs_truthLead_" + cls + "_pTgammaTruth_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];#DeltaR",
                                cfg.ptBins, 70, 0.0, 3.5))
              h->Fill(leadTruth->pt, std::isfinite(truthLeadDR) ? truthLeadDR : dR(recoLead->eta, recoLead->phi, truthLead->eta, truthLead->phi), ev.weight);
            if (auto* h = hist2Uniform(out, trig, "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_" + cls + "_" + rKey + cs,
                                       "h2_leadTruthRecoilMatch_xJRecoJet1_vs_dphiRecoJet1_" + cls + "_" + rKey + cs + ";|#Delta#phi|;x_{J}^{reco}",
                                       64, 0.0, M_PI, 60, 0.0, 3.0))
              h->Fill(dphiReco, xJReco, ev.weight);
            if (truthLeadMatchedReco)
            {
              if (auto* h = hist2Uniform(out, trig, "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_" + cls + "_" + rKey + cs,
                                         "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTtruthLead_" + cls + "_" + rKey + cs + ";p_{T}^{lead truth recoil} [GeV];p_{T}^{recoilJet1,reco} [GeV]",
                                         120, 0.0, 60.0, 120, 0.0, 60.0))
                h->Fill(truthLead->pt, recoLead->pt, ev.weight);
              if (auto* h = hist2Uniform(out, trig, "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_" + cls + "_" + rKey + cs,
                                         "h2_leadTruthRecoilMatch_pTrecoJet1_vs_pTrecoTruthMatch_" + cls + "_" + rKey + cs + ";p_{T}^{reco match} [GeV];p_{T}^{recoilJet1,reco} [GeV]",
                                         120, 0.0, 60.0, 120, 0.0, 60.0))
                h->Fill(truthLeadMatchedReco->pt, recoLead->pt, ev.weight);
              if (auto* h = hist2Uniform(out, trig, "h2_leadRecoilJetPtResponse_pTtruth_ratio_" + rKey + cs,
                                         "h2_leadRecoilJetPtResponse_pTtruth_ratio_" + rKey + cs + ";p_{T}^{jet,truth} [GeV];p_{T}^{jet,reco}/p_{T}^{jet,truth}",
                                         120, 0.0, 60.0, 120, 0.0, 2.0))
                h->Fill(truthLead->pt, truthLeadMatchedReco->pt / truthLead->pt, ev.weight);
            }
          }

          const double xJTruth = truthLead->pt / leadTruth->pt;
          double subTruthPt = 0.0;
          for (const JetRec* tj : truthRecoil)
            if (tj != truthLead && tj->pt > subTruthPt) subTruthPt = tj->pt;
          const double alphaTruth = subTruthPt > 0.0 ? subTruthPt / leadTruth->pt : 0.0;
          if (leadReco && leadReco->truthSignal)
          {
            if (num)
            {
              if (auto* h = hist3VarUniform(out, trig, "h_JES3Truth_pT_xJ_alpha_" + rKey + cs,
                                            "h_JES3Truth_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J}^{truth};#alpha^{truth}",
                                            cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
                h->Fill(leadTruth->pt, xJTruth, alphaTruth, ev.weight);
              if (auto* h = hist3VarUniform(out, trig, "h_JES3Truth_pT_jet1Pt_alpha_" + rKey + cs,
                                            "h_JES3Truth_pT_jet1Pt_alpha_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];p_{T}^{jet1,truth} [GeV];#alpha^{truth}",
                                            cfg.ptBins, 120, 0.0, 60.0, 40, 0.0, 2.0))
                h->Fill(leadTruth->pt, truthLead->pt, alphaTruth, ev.weight);
              if (recoLead)
              {
                const double recoAlpha = 0.0;
                if (auto* h = hist3VarUniform(out, trig, "h_JES3RecoTruthTagged_pT_xJ_alpha_" + rKey + cs,
                                              "h_JES3RecoTruthTagged_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma,reco} [GeV];x_{J}^{reco};#alpha^{reco}",
                                              cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
                  h->Fill(leadReco->pt, recoLead->pt / leadReco->pt, recoAlpha, ev.weight);
              }
            }
            else
            {
              if (auto* h = hist3VarUniform(out, trig, "h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_" + rKey + cs,
                                            "h_JES3TruthRecoCondNoJetMatch_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma,truth} [GeV];x_{J}^{truth};#alpha^{truth}",
                                            cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
                h->Fill(leadTruth->pt, xJTruth, alphaTruth, ev.weight);
            }
          }
        }
      }
    }
  }
}

void Fun4All_recoilJets_poolReplay(const int nEvents,
                                   const char* fileList,
                                   const char* outRoot,
                                   bool /*unused*/)
{
  using namespace rjpool;

  Config cfg = loadConfig();
  std::vector<Entry> entries = gForcedFanoutEntries ? *gForcedFanoutEntries : loadFanout(outRoot ? outRoot : "poolReplay.root", cfg);
  const bool exportTrainingMode = envBool("RJ_POOL_EXPORT_AUAU_BDT_TRAINING_TREE", false);
  const bool strictLegacyOutput = envBool("RJ_REPLAY_STRICT_LEGACY_OUTPUT", true);
  const bool writeViewCatalog = envBool("RJ_REPLAY_WRITE_VIEW_CATALOG", false);
  TH1::AddDirectory(false);
  const std::size_t maxOpenOutputs = replayMaxOpenOutputs();
  std::vector<std::string> uniqueRoots;
  std::unordered_map<std::string, int> rootSeen;
  for (const Entry& e : entries)
  {
    if (!rootSeen.count(e.outRoot))
    {
      rootSeen[e.outRoot] = 1;
      uniqueRoots.push_back(e.outRoot);
    }
  }
  if (!strictLegacyOutput && !gForcedFanoutEntries && uniqueRoots.size() > maxOpenOutputs)
  {
    std::cerr << "[RecoilJetsPoolReplay] Fanout has " << uniqueRoots.size()
              << " outputs; batching with RJ_REPLAY_MAX_OPEN_OUTPUTS="
              << maxOpenOutputs << "\n";
    for (std::size_t first = 0; first < uniqueRoots.size(); first += maxOpenOutputs)
    {
      const std::size_t last = std::min(uniqueRoots.size(), first + maxOpenOutputs);
      std::unordered_map<std::string, int> batchRoots;
      for (std::size_t i = first; i < last; ++i) batchRoots[uniqueRoots[i]] = 1;
      std::vector<Entry> batch;
      for (const Entry& e : entries)
      {
        if (batchRoots.count(e.outRoot)) batch.push_back(e);
      }
      gForcedFanoutEntries = &batch;
      Fun4All_recoilJets_poolReplay(nEvents, fileList, outRoot, false);
      gForcedFanoutEntries = nullptr;
    }
    return;
  }

  TChain events("AnalysisEventPool");
  TChain photons("AnalysisPhotonPool");
  TChain jets("AnalysisJetPool");
  TChain truthPhotons("AnalysisTruthPhotonPool");

  std::ifstream in(fileList);
  if (in)
  {
    for (std::string line; std::getline(in, line); )
    {
      line = trim(line);
      if (line.empty() || line[0] == '#') continue;
      events.Add(line.c_str());
      photons.Add(line.c_str());
      jets.Add(line.c_str());
      truthPhotons.Add(line.c_str());
    }
  }
  else
  {
    events.Add(fileList);
    photons.Add(fileList);
    jets.Add(fileList);
    truthPhotons.Add(fileList);
  }

  if (events.GetEntries() <= 0 || photons.GetEntries() < 0)
  {
    std::cerr << "[RecoilJetsPoolReplay] No pool entries found from " << fileList << "\n";
    return;
  }

  const bool schemaOk =
    validateBranches(events, "AnalysisEventPool",
                     {"schema", "run", "evt", "eventKey", "isSim", "isAuAu", "centBin", "centPercent", "vz", "weight", "triggers"}, true) &&
    validateBranches(photons, "AnalysisPhotonPool",
                     {"schema", "eventKey", "pt", "eta", "phi", "energy", "eiso", "eiso_r03", "eiso_r04",
                      "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx", "weta35_cogx", "wphi53_cogx",
                      "et1", "et2", "et3", "et4", "e11_over_e33", "e32_over_e35",
                      "iso_03_emcal", "iso_03_hcalin", "iso_03_hcalout",
                      "iso_04_emcal", "iso_04_hcalin", "iso_04_hcalout",
                      "e11_over_e22", "e11_over_e13", "e11_over_e15", "e11_over_e17", "e11_over_e31",
                      "e11_over_e51", "e11_over_e71", "e22_over_e33", "e22_over_e35", "e22_over_e37",
                      "e22_over_e53", "w32", "w52", "w72", "mean_time", "npb_score", "tight_bdt_score",
                      "auau_npb_score", "auau_tight_bdt_score", "mbd_time", "cluster_mbd_delta_t",
                      "npb_has_away_jet", "npb_label", "is_npb", "truthSignal", "truthTrackId",
                      "truthPt", "truthEta", "truthPhi"}, true) &&
    validateBranches(jets, "AnalysisJetPool",
                     {"schema", "eventKey", "rKey", "isTruth", "pt", "raw_pt", "areaSub_pt", "eta", "phi", "jet_area"}, false) &&
    validateBranches(truthPhotons, "AnalysisTruthPhotonPool",
                     {"schema", "eventKey", "trackId", "barcode", "pt", "eta", "phi", "iso"}, false);
  if (!schemaOk) return;

  std::vector<Output> outputs;
  std::unordered_map<std::string, std::size_t> outputIndex;
  if (!exportTrainingMode) for (const Entry& e : entries)
  {
    const auto found = outputIndex.find(e.outRoot);
    std::size_t idx = 0;
    if (found == outputIndex.end())
    {
      Output o;
      o.outRoot = e.outRoot;
      outputs.push_back(o);
      idx = outputs.size() - 1;
      outputIndex[e.outRoot] = idx;
    }
    else
    {
      idx = found->second;
    }
    View v;
    v.entry = e;
    v.cfg = e.cfg;
    v.viewId = e.viewId.empty() ? defaultViewId(e.cfg) : e.viewId;
    v.materialize = e.materialize.empty() ? "physics" : e.materialize;
    outputs[idx].views.push_back(v);
  }
  for (Output& o : outputs)
  {
    const bool singleView = (o.views.size() == 1);
    for (View& v : o.views)
    {
      v.legacyTopLevel = strictLegacyOutput || singleView || v.viewId == "legacy" || v.viewId == "nominal_legacy";
    }
  }

  int eventRun = 0;
  int eventEvt = 0;
  long long eventKey = 0;
  int isSim = 0;
  int isAuAu = 0;
  int centBin = -1;
  float centPercent = -1.0f;
  float vz = 0.0f;
  float weight = 1.0f;
  std::vector<std::string>* triggers = nullptr;
  events.SetBranchAddress("run", &eventRun);
  events.SetBranchAddress("evt", &eventEvt);
  events.SetBranchAddress("eventKey", &eventKey);
  events.SetBranchAddress("isSim", &isSim);
  events.SetBranchAddress("isAuAu", &isAuAu);
  events.SetBranchAddress("centBin", &centBin);
  events.SetBranchAddress("centPercent", &centPercent);
  events.SetBranchAddress("vz", &vz);
  events.SetBranchAddress("weight", &weight);
  events.SetBranchAddress("triggers", &triggers);

  std::unordered_map<long long, EventRec> eventMap;
  std::vector<long long> eventOrder;
  const Long64_t nEvt = events.GetEntries();
  for (Long64_t i = 0; i < nEvt; ++i)
  {
    events.GetEntry(i);
    if (nEvents > 0 && static_cast<int>(i) >= nEvents) break;
    const long long key = scopedEventKey(events, eventKey);
    EventRec rec;
    rec.run = eventRun;
    rec.evt = eventEvt;
    rec.isSim = (isSim != 0);
    rec.isAuAu = (isAuAu != 0);
    rec.centBin = centBin;
    rec.centPercent = centPercent;
    rec.vz = vz;
    rec.weight = weight;
    if (triggers) rec.triggers = *triggers;
    if (rec.triggers.empty()) rec.triggers.push_back(rec.isAuAu ? "ALL" : "SIM");
    eventMap[key] = rec;
    eventOrder.push_back(key);

    for (Output& o : outputs)
    {
      for (View& view : o.views)
      {
        o.activeView = &view;
        const Config& ocfg = view.cfg;
        if (ocfg.useVzCut && std::fabs(static_cast<double>(rec.vz)) >= ocfg.vzCutCm) continue;
        for (const std::string& trig : rec.triggers)
        {
          countHist(o, trig, "cnt_" + trig)->Fill(1, rec.weight);
          vertexHist(o, trig)->Fill(rec.vz, rec.weight);
        }
      }
      o.activeView = nullptr;
    }
  }

  std::unordered_map<long long, std::vector<PhotonRec>> photonsByEvent;
  float pt = 0, eta = 0, phi = 0, energy = 0, eiso = 0, eisoR03 = 0, eisoR04 = 0;
  float iso03Emcal = 0, iso03Hcalin = 0, iso03Hcalout = 0, iso04Emcal = 0, iso04Hcalin = 0, iso04Hcalout = 0;
  float weta = 0, wphi = 0, weta33 = 0, wphi33 = 0, weta35 = 0, wphi53 = 0;
  float et1 = 0, et2 = 0, et3 = 0, et4 = 0, e11e33 = 0, e32e35 = 0;
  float e11e22 = 0, e11e13 = 0, e11e15 = 0, e11e17 = 0, e11e31 = 0, e11e51 = 0, e11e71 = 0;
  float e22e33 = 0, e22e35 = 0, e22e37 = 0, e22e53 = 0, w32 = 0, w52 = 0, w72 = 0, meanTime = 0;
  float npb = 0, tightBDT = 0, auauNpb = 0, auauTightBDT = 0;
  float mbdTime = 0, clusterMbdDeltaT = 0;
  int phoIndex = -1, npbHasAwayJet = 0, npbLabel = -1, isNPB = -1;
  int truthSignal = 0, truthTrackId = -1;
  float truthPt = 0, truthEta = 0, truthPhi = 0;
  std::vector<float>* extraFeatureValues = nullptr;
  photons.SetBranchAddress("eventKey", &eventKey);
  if (photons.GetBranch("index")) photons.SetBranchAddress("index", &phoIndex);
  photons.SetBranchAddress("pt", &pt);
  photons.SetBranchAddress("eta", &eta);
  photons.SetBranchAddress("phi", &phi);
  photons.SetBranchAddress("energy", &energy);
  photons.SetBranchAddress("eiso", &eiso);
  photons.SetBranchAddress("eiso_r03", &eisoR03);
  photons.SetBranchAddress("eiso_r04", &eisoR04);
  photons.SetBranchAddress("iso_03_emcal", &iso03Emcal);
  photons.SetBranchAddress("iso_03_hcalin", &iso03Hcalin);
  photons.SetBranchAddress("iso_03_hcalout", &iso03Hcalout);
  photons.SetBranchAddress("iso_04_emcal", &iso04Emcal);
  photons.SetBranchAddress("iso_04_hcalin", &iso04Hcalin);
  photons.SetBranchAddress("iso_04_hcalout", &iso04Hcalout);
  photons.SetBranchAddress("weta_cogx", &weta);
  photons.SetBranchAddress("wphi_cogx", &wphi);
  photons.SetBranchAddress("weta33_cogx", &weta33);
  photons.SetBranchAddress("wphi33_cogx", &wphi33);
  photons.SetBranchAddress("weta35_cogx", &weta35);
  photons.SetBranchAddress("wphi53_cogx", &wphi53);
  photons.SetBranchAddress("et1", &et1);
  photons.SetBranchAddress("et2", &et2);
  photons.SetBranchAddress("et3", &et3);
  photons.SetBranchAddress("et4", &et4);
  photons.SetBranchAddress("e11_over_e33", &e11e33);
  photons.SetBranchAddress("e32_over_e35", &e32e35);
  photons.SetBranchAddress("e11_over_e22", &e11e22);
  photons.SetBranchAddress("e11_over_e13", &e11e13);
  photons.SetBranchAddress("e11_over_e15", &e11e15);
  photons.SetBranchAddress("e11_over_e17", &e11e17);
  photons.SetBranchAddress("e11_over_e31", &e11e31);
  photons.SetBranchAddress("e11_over_e51", &e11e51);
  photons.SetBranchAddress("e11_over_e71", &e11e71);
  photons.SetBranchAddress("e22_over_e33", &e22e33);
  photons.SetBranchAddress("e22_over_e35", &e22e35);
  photons.SetBranchAddress("e22_over_e37", &e22e37);
  photons.SetBranchAddress("e22_over_e53", &e22e53);
  photons.SetBranchAddress("w32", &w32);
  photons.SetBranchAddress("w52", &w52);
  photons.SetBranchAddress("w72", &w72);
  photons.SetBranchAddress("mean_time", &meanTime);
  photons.SetBranchAddress("npb_score", &npb);
  photons.SetBranchAddress("tight_bdt_score", &tightBDT);
  photons.SetBranchAddress("auau_npb_score", &auauNpb);
  photons.SetBranchAddress("auau_tight_bdt_score", &auauTightBDT);
  photons.SetBranchAddress("mbd_time", &mbdTime);
  photons.SetBranchAddress("cluster_mbd_delta_t", &clusterMbdDeltaT);
  photons.SetBranchAddress("npb_has_away_jet", &npbHasAwayJet);
  photons.SetBranchAddress("npb_label", &npbLabel);
  photons.SetBranchAddress("is_npb", &isNPB);
  photons.SetBranchAddress("truthSignal", &truthSignal);
  photons.SetBranchAddress("truthTrackId", &truthTrackId);
  photons.SetBranchAddress("truthPt", &truthPt);
  photons.SetBranchAddress("truthEta", &truthEta);
  photons.SetBranchAddress("truthPhi", &truthPhi);
  if (photons.GetBranch("extra_feature_values"))
  {
    photons.SetBranchAddress("extra_feature_values", &extraFeatureValues);
  }

  const Long64_t nPho = photons.GetEntries();
  for (Long64_t i = 0; i < nPho; ++i)
  {
    photons.GetEntry(i);
    const long long key = scopedEventKey(photons, eventKey);
    const auto evIt = eventMap.find(key);
    if (evIt == eventMap.end()) continue;
    PhotonRec pr;
    pr.eventKey = key;
    pr.index = phoIndex;
    pr.pt = pt;
    pr.eta = eta;
    pr.phi = phi;
    pr.energy = energy;
    pr.eiso = eiso;
    pr.eisoR03 = eisoR03;
    pr.eisoR04 = eisoR04;
    pr.iso03Emcal = iso03Emcal;
    pr.iso03Hcalin = iso03Hcalin;
    pr.iso03Hcalout = iso03Hcalout;
    pr.iso04Emcal = iso04Emcal;
    pr.iso04Hcalin = iso04Hcalin;
    pr.iso04Hcalout = iso04Hcalout;
    pr.weta = weta;
    pr.wphi = wphi;
    pr.weta33 = weta33;
    pr.wphi33 = wphi33;
    pr.weta35 = weta35;
    pr.wphi53 = wphi53;
    pr.et1 = et1;
    pr.et2 = et2;
    pr.et3 = et3;
    pr.et4 = et4;
    pr.e11e33 = e11e33;
    pr.e32e35 = e32e35;
    pr.e11e22 = e11e22;
    pr.e11e13 = e11e13;
    pr.e11e15 = e11e15;
    pr.e11e17 = e11e17;
    pr.e11e31 = e11e31;
    pr.e11e51 = e11e51;
    pr.e11e71 = e11e71;
    pr.e22e33 = e22e33;
    pr.e22e35 = e22e35;
    pr.e22e37 = e22e37;
    pr.e22e53 = e22e53;
    pr.w32 = w32;
    pr.w52 = w52;
    pr.w72 = w72;
    pr.meanTime = meanTime;
    pr.npb = npb;
    pr.tightBDT = tightBDT;
    pr.auauNpb = auauNpb;
    pr.auauTightBDT = auauTightBDT;
    pr.mbdTime = mbdTime;
    pr.clusterMbdDeltaT = clusterMbdDeltaT;
    pr.npbHasAwayJet = npbHasAwayJet;
    pr.npbLabel = npbLabel;
    pr.isNPB = isNPB;
    pr.truthSignal = truthSignal;
    pr.truthTrackId = truthTrackId;
    pr.truthPt = truthPt;
    pr.truthEta = truthEta;
    pr.truthPhi = truthPhi;
    if (extraFeatureValues) pr.extraFeatures = *extraFeatureValues;
    photonsByEvent[key].push_back(pr);
  }

  std::unordered_map<long long, std::vector<JetRec>> jetsByEvent;
  long long jetEventKey = 0;
  int jetIsTruth = 0;
  float jetPt = 0, jetRawPt = 0, jetAreaSubPt = 0, jetEta = 0, jetPhi = 0, jetMass = 0, jetArea = 0;
  std::string* jetRKey = nullptr;
  if (jets.GetEntries() > 0)
  {
    jets.SetBranchAddress("eventKey", &jetEventKey);
    jets.SetBranchAddress("rKey", &jetRKey);
    jets.SetBranchAddress("isTruth", &jetIsTruth);
    jets.SetBranchAddress("pt", &jetPt);
    jets.SetBranchAddress("raw_pt", &jetRawPt);
    jets.SetBranchAddress("areaSub_pt", &jetAreaSubPt);
    jets.SetBranchAddress("eta", &jetEta);
    jets.SetBranchAddress("phi", &jetPhi);
    if (jets.GetBranch("mass")) jets.SetBranchAddress("mass", &jetMass);
    jets.SetBranchAddress("jet_area", &jetArea);
    const Long64_t nJet = jets.GetEntries();
    for (Long64_t i = 0; i < nJet; ++i)
    {
      jets.GetEntry(i);
      const long long key = scopedEventKey(jets, jetEventKey);
      if (eventMap.find(key) == eventMap.end()) continue;
      JetRec jr;
      jr.eventKey = key;
      jr.rKey = jetRKey ? *jetRKey : "";
      jr.isTruth = jetIsTruth;
      jr.pt = jetPt;
      jr.rawPt = jetRawPt;
      jr.areaSubPt = jetAreaSubPt;
      jr.eta = jetEta;
      jr.phi = jetPhi;
      jr.mass = jets.GetBranch("mass") ? jetMass : std::numeric_limits<float>::quiet_NaN();
      jr.area = jetArea;
      if (!jr.rKey.empty()) jetsByEvent[key].push_back(jr);
    }
  }

  std::unordered_map<long long, std::vector<TruthPhotonRec>> truthByEvent;
  long long truthEventKey = 0;
  int truthIdx = 0, truthBarcode = -1;
  float tpt = 0, teta = 0, tphi = 0, tiso = 0;
  if (truthPhotons.GetEntries() > 0)
  {
    truthPhotons.SetBranchAddress("eventKey", &truthEventKey);
    truthPhotons.SetBranchAddress("trackId", &truthIdx);
    truthPhotons.SetBranchAddress("barcode", &truthBarcode);
    truthPhotons.SetBranchAddress("pt", &tpt);
    truthPhotons.SetBranchAddress("eta", &teta);
    truthPhotons.SetBranchAddress("phi", &tphi);
    truthPhotons.SetBranchAddress("iso", &tiso);
    const Long64_t nTruthPho = truthPhotons.GetEntries();
    for (Long64_t i = 0; i < nTruthPho; ++i)
    {
      truthPhotons.GetEntry(i);
      const long long key = scopedEventKey(truthPhotons, truthEventKey);
      if (eventMap.find(key) == eventMap.end()) continue;
      TruthPhotonRec tr;
      tr.eventKey = key;
      tr.trackId = truthIdx;
      tr.barcode = truthBarcode;
      tr.pt = tpt;
      tr.eta = teta;
      tr.phi = tphi;
      tr.iso = tiso;
      if (std::isfinite(tr.pt) && tr.pt > 0.0) truthByEvent[key].push_back(tr);
    }
  }

  if (exportTrainingMode)
  {
    const char* trainOutEnv = std::getenv("RJ_POOL_AUAU_BDT_TRAINING_OUT");
    const std::string trainOut =
      (trainOutEnv && *trainOutEnv) ? std::string(trainOutEnv)
                                    : std::string(outRoot && *outRoot ? outRoot : "AuAuPhotonIDTrainingFromPool.root");
    TFile trainFile(trainOut.c_str(), "RECREATE");
    if (trainFile.IsZombie())
    {
      std::cerr << "[RecoilJetsPoolReplay] Could not open training export output " << trainOut << "\n";
      return;
    }

    int run = 0;
    long long evt = 0;
    int is_signal = 0;
    int pt_bin = -1;
    int cent_bin = -1;
    float cluster_Et = 0, cluster_Eta = 0, cluster_Phi = 0, centrality = -1, vertexz = 0, event_weight = 1;
    float reco_eiso = 0, cluster_mean_time = -999, mbd_time_out = -999, cluster_mbd_delta_t = -999;
    int npb_label = -1, is_npb = -1, npb_has_away_jet = 0;
    float cluster_weta_cogx = 0, cluster_wphi_cogx = 0, cluster_weta33_cogx = 0, cluster_wphi33_cogx = 0;
    float cluster_weta35_cogx = 0, cluster_wphi53_cogx = 0;
    float cluster_et1 = 0, cluster_et2 = 0, cluster_et3 = 0, cluster_et4 = 0;
    float e11_over_e33 = 0, e32_over_e35 = 0, e11_over_e22 = 0, e11_over_e13 = 0, e11_over_e15 = 0;
    float e11_over_e17 = 0, e11_over_e31 = 0, e11_over_e51 = 0, e11_over_e71 = 0;
    float e22_over_e33 = 0, e22_over_e35 = 0, e22_over_e37 = 0, e22_over_e53 = 0;
    float cluster_w32 = 0, cluster_w52 = 0, cluster_w72 = 0;
    float npb_score = -2, auau_npb_score = -2, auau_tight_bdt_score = -2;

    TTree tree("AuAuPhotonIDTrainingTree", "AuAu photon-ID/NPB BDT training candidates exported from AnalysisPool");
    tree.Branch("run", &run, "run/I");
    tree.Branch("evt", &evt, "evt/L");
    tree.Branch("is_signal", &is_signal, "is_signal/I");
    tree.Branch("pt_bin", &pt_bin, "pt_bin/I");
    tree.Branch("cent_bin", &cent_bin, "cent_bin/I");
    tree.Branch("cluster_Et", &cluster_Et, "cluster_Et/F");
    tree.Branch("cluster_Eta", &cluster_Eta, "cluster_Eta/F");
    tree.Branch("cluster_Phi", &cluster_Phi, "cluster_Phi/F");
    tree.Branch("centrality", &centrality, "centrality/F");
    tree.Branch("vertexz", &vertexz, "vertexz/F");
    tree.Branch("event_weight", &event_weight, "event_weight/F");
    tree.Branch("reco_eiso", &reco_eiso, "reco_eiso/F");
    tree.Branch("npb_label", &npb_label, "npb_label/I");
    tree.Branch("is_npb", &is_npb, "is_npb/I");
    tree.Branch("cluster_mean_time", &cluster_mean_time, "cluster_mean_time/F");
    tree.Branch("mbd_time", &mbd_time_out, "mbd_time/F");
    tree.Branch("cluster_mbd_delta_t", &cluster_mbd_delta_t, "cluster_mbd_delta_t/F");
    tree.Branch("npb_has_away_jet", &npb_has_away_jet, "npb_has_away_jet/I");
    tree.Branch("cluster_weta_cogx", &cluster_weta_cogx, "cluster_weta_cogx/F");
    tree.Branch("cluster_wphi_cogx", &cluster_wphi_cogx, "cluster_wphi_cogx/F");
    tree.Branch("cluster_weta33_cogx", &cluster_weta33_cogx, "cluster_weta33_cogx/F");
    tree.Branch("cluster_wphi33_cogx", &cluster_wphi33_cogx, "cluster_wphi33_cogx/F");
    tree.Branch("cluster_weta35_cogx", &cluster_weta35_cogx, "cluster_weta35_cogx/F");
    tree.Branch("cluster_wphi53_cogx", &cluster_wphi53_cogx, "cluster_wphi53_cogx/F");
    tree.Branch("cluster_et1", &cluster_et1, "cluster_et1/F");
    tree.Branch("cluster_et2", &cluster_et2, "cluster_et2/F");
    tree.Branch("cluster_et3", &cluster_et3, "cluster_et3/F");
    tree.Branch("cluster_et4", &cluster_et4, "cluster_et4/F");
    tree.Branch("e11_over_e33", &e11_over_e33, "e11_over_e33/F");
    tree.Branch("e32_over_e35", &e32_over_e35, "e32_over_e35/F");
    tree.Branch("e11_over_e22", &e11_over_e22, "e11_over_e22/F");
    tree.Branch("e11_over_e13", &e11_over_e13, "e11_over_e13/F");
    tree.Branch("e11_over_e15", &e11_over_e15, "e11_over_e15/F");
    tree.Branch("e11_over_e17", &e11_over_e17, "e11_over_e17/F");
    tree.Branch("e11_over_e31", &e11_over_e31, "e11_over_e31/F");
    tree.Branch("e11_over_e51", &e11_over_e51, "e11_over_e51/F");
    tree.Branch("e11_over_e71", &e11_over_e71, "e11_over_e71/F");
    tree.Branch("e22_over_e33", &e22_over_e33, "e22_over_e33/F");
    tree.Branch("e22_over_e35", &e22_over_e35, "e22_over_e35/F");
    tree.Branch("e22_over_e37", &e22_over_e37, "e22_over_e37/F");
    tree.Branch("e22_over_e53", &e22_over_e53, "e22_over_e53/F");
    tree.Branch("cluster_w32", &cluster_w32, "cluster_w32/F");
    tree.Branch("cluster_w52", &cluster_w52, "cluster_w52/F");
    tree.Branch("cluster_w72", &cluster_w72, "cluster_w72/F");
    tree.Branch("npb_score", &npb_score, "npb_score/F");
    tree.Branch("auau_npb_score", &auau_npb_score, "auau_npb_score/F");
    tree.Branch("auau_tight_bdt_score", &auau_tight_bdt_score, "auau_tight_bdt_score/F");

    const Config& trainCfg = entries.empty() ? cfg : entries.front().cfg;
    for (long long key : eventOrder)
    {
      const auto evIt = eventMap.find(key);
      if (evIt == eventMap.end()) continue;
      const EventRec& ev = evIt->second;
      const std::vector<PhotonRec>& evPhotons = photonsByEvent[key];
      for (const PhotonRec& ph : evPhotons)
      {
        if (!std::isfinite(ph.pt) || ph.pt <= 0.0) continue;
        run = ev.run;
        evt = ev.evt;
        is_signal = ph.truthSignal ? 1 : 0;
        pt_bin = findPtBin(ph.pt, trainCfg.ptBins);
        cent_bin = ev.isAuAu ? findCentBin(ev.centBin, trainCfg.centEdges) : -1;
        cluster_Et = ph.pt;
        cluster_Eta = ph.eta;
        cluster_Phi = ph.phi;
        centrality = std::isfinite(ev.centPercent) ? ev.centPercent : static_cast<float>(ev.centBin);
        vertexz = ev.vz;
        event_weight = ev.weight;
        reco_eiso = std::isfinite(ph.eisoR04) ? ph.eisoR04 : ph.eiso;
        npb_label = ph.npbLabel;
        is_npb = ph.isNPB;
        cluster_mean_time = std::isfinite(ph.meanTime) ? ph.meanTime : -999.0f;
        mbd_time_out = std::isfinite(ph.mbdTime) ? ph.mbdTime : -999.0f;
        cluster_mbd_delta_t = std::isfinite(ph.clusterMbdDeltaT) ? ph.clusterMbdDeltaT : -999.0f;
        npb_has_away_jet = ph.npbHasAwayJet;
        cluster_weta_cogx = ph.weta;
        cluster_wphi_cogx = ph.wphi;
        cluster_weta33_cogx = ph.weta33;
        cluster_wphi33_cogx = ph.wphi33;
        cluster_weta35_cogx = ph.weta35;
        cluster_wphi53_cogx = ph.wphi53;
        cluster_et1 = ph.et1;
        cluster_et2 = ph.et2;
        cluster_et3 = ph.et3;
        cluster_et4 = ph.et4;
        e11_over_e33 = ph.e11e33;
        e32_over_e35 = ph.e32e35;
        e11_over_e22 = ph.e11e22;
        e11_over_e13 = ph.e11e13;
        e11_over_e15 = ph.e11e15;
        e11_over_e17 = ph.e11e17;
        e11_over_e31 = ph.e11e31;
        e11_over_e51 = ph.e11e51;
        e11_over_e71 = ph.e11e71;
        e22_over_e33 = ph.e22e33;
        e22_over_e35 = ph.e22e35;
        e22_over_e37 = ph.e22e37;
        e22_over_e53 = ph.e22e53;
        cluster_w32 = ph.w32;
        cluster_w52 = ph.w52;
        cluster_w72 = ph.w72;
        npb_score = ph.npb;
        auau_npb_score = ph.auauNpb;
        auau_tight_bdt_score = ph.auauTightBDT;
        tree.Fill();
      }
    }
    tree.Write("", TObject::kOverwrite);
    trainFile.Close();
    std::cout << "[RecoilJetsPoolReplay] Exported AuAuPhotonIDTrainingTree to "
              << trainOut << " entries=" << tree.GetEntries() << "\n";
    return;
  }

  ReplayBDTModels replayBdt;

  for (long long key : eventOrder)
  {
    const auto evIt = eventMap.find(key);
    if (evIt == eventMap.end()) continue;
    const EventRec& ev = evIt->second;

    const std::vector<PhotonRec>& evPhotons = photonsByEvent[key];
    const std::vector<JetRec>& evJets = jetsByEvent[key];
    const std::vector<TruthPhotonRec>& evTruth = truthByEvent[key];
    const TruthPhotonRec* leadTruth = nullptr;
    for (const TruthPhotonRec& tp : evTruth)
    {
      if (!leadTruth || tp.pt > leadTruth->pt) leadTruth = &tp;
    }

    for (Output& o : outputs)
    {
      for (View& view : o.views)
      {
        o.activeView = &view;
        const Config& ocfg = view.cfg;
        if (ocfg.useVzCut && std::fabs(static_cast<double>(ev.vz)) >= ocfg.vzCutCm) continue;
        const int centIdx = ev.isAuAu ? findCentBin(ev.centBin, ocfg.centEdges) : -1;
        if (ev.isAuAu && centIdx < 0) continue;

        fillInclusiveJetQA(o, ev, centIdx, ocfg, evJets);

        const PhotonRec* leadIsoTight = nullptr;
        const PhotonRec* leadIsoNonTight = nullptr;
        const PhotonRec* leadNonIsoTight = nullptr;
        const PhotonRec* leadNonIsoNonTight = nullptr;

        for (const PhotonRec& ph : evPhotons)
        {
          if (!std::isfinite(ph.pt) || ph.pt <= 0.0) continue;
          if (!std::isfinite(ph.eta) || std::fabs(ph.eta) >= 0.7) continue;
          const int ptIdx = findPtBin(ph.pt, ocfg.ptBins);
          if (ptIdx < 0) continue;

          const double isoValue = photonIsoForView(ph, ocfg);
          const double isoThr = ocfg.slidingIso ? (ocfg.isoA + ocfg.isoB * ph.pt) : ocfg.isoFixed;
          const double nonIsoThr = isoThr + ocfg.isoGap;
          const bool iso = std::isfinite(isoValue) && isoValue < isoThr;
          const bool nonIso = std::isfinite(isoValue) && isoValue > nonIsoThr;

          const double auauNpbScore = replayAuAuNpbScore(replayBdt, view.entry, ev, ph);
          const double auauTightBDTScore = replayAuAuTightScore(replayBdt, view.entry, ev, ph);
          if (!passPre(view.entry, ocfg, ph.pt, ph.weta, ph.et1, ph.e11e33, ph.e32e35, ph.npb, auauNpbScore)) continue;
          const int tag = tightTag(view.entry, ocfg, ph.pt, ph.weta, ph.wphi, ph.et1, ph.e11e33, ph.e32e35, ph.tightBDT, auauTightBDTScore);

          const char* base = nullptr;
          int regionBin = 0;
          if (iso && tag == 1) { base = "h_isIsolated_isTight"; regionBin = 1; }
          else if (nonIso && tag == 1) { base = "h_notIsolated_isTight"; regionBin = 2; }
          else if (iso && tag == 2) { base = "h_isIsolated_notTight"; regionBin = 3; }
          else if (nonIso && tag == 2) { base = "h_notIsolated_notTight"; regionBin = 4; }

          const std::string sfx = suffix(ptIdx, centIdx, ev.isAuAu, ocfg);
          fillPhotonLegacyQA(o, ev, centIdx, ocfg, ph, tag, iso, nonIso, regionBin);

          if (tag == 3 || !base) continue;

          for (const std::string& trig : ev.triggers)
          {
            if (auto* h = abcdHist(o, trig, base, sfx)) h->Fill(1.0, ev.weight);
            if (ph.truthSignal && regionBin > 0)
            {
              if (auto* hSig = sigABCDHist(o, trig, sfx)) hSig->Fill(regionBin, ev.weight);
            }
          }

          if (iso && tag == 1)
          {
            if (!leadIsoTight || ph.pt > leadIsoTight->pt) leadIsoTight = &ph;
          }
          else if (iso && tag == 2)
          {
            if (!leadIsoNonTight || ph.pt > leadIsoNonTight->pt) leadIsoNonTight = &ph;
          }
          else if (nonIso && tag == 1)
          {
            if (!leadNonIsoTight || ph.pt > leadNonIsoTight->pt) leadNonIsoTight = &ph;
          }
          else if (nonIso && tag == 2)
          {
            if (!leadNonIsoNonTight || ph.pt > leadNonIsoNonTight->pt) leadNonIsoNonTight = &ph;
          }
        }

        auto fillLeadPurity = [&](const PhotonRec* ph, const std::string& base, int regionBin)
        {
          if (!ph) return;
          const int pidx = findPtBin(ph->pt, ocfg.ptBins);
          if (pidx < 0) return;
          const std::string sfx = suffix(pidx, centIdx, ev.isAuAu, ocfg);
          for (const std::string& trig : ev.triggers)
          {
            if (auto* h = abcdHist(o, trig, "h_xJpurityLead_" + base, sfx)) h->Fill(1.0, ev.weight);
            if (ph->truthSignal && regionBin > 0)
            {
              const std::string name = "h_xJpurityLead_sigABCD_MC" + sfx;
              auto& m = o.h[histScopeKey(o, trig)];
              TH1I* hSig = nullptr;
              if (auto it = m.find(name); it != m.end()) hSig = dynamic_cast<TH1I*>(it->second);
              if (!hSig)
              {
                ensureDir(o, trig);
                hSig = new TH1I(name.c_str(), (name + ";Leading reco ABCD region for matched truth-signal photons;Entries").c_str(), 4, 0.5, 4.5);
                hSig->GetXaxis()->SetBinLabel(1, "A");
                hSig->GetXaxis()->SetBinLabel(2, "B");
                hSig->GetXaxis()->SetBinLabel(3, "C");
                hSig->GetXaxis()->SetBinLabel(4, "D");
                hSig->Sumw2();
                m[name] = hSig;
              }
              hSig->Fill(regionBin, ev.weight);
            }
          }
        };
        fillLeadPurity(leadIsoTight, "isIsolated_isTight", 1);
        fillLeadPurity(leadNonIsoTight, "notIsolated_isTight", 2);
        fillLeadPurity(leadIsoNonTight, "isIsolated_notTight", 3);
        fillLeadPurity(leadNonIsoNonTight, "notIsolated_notTight", 4);

        fillPhotonCounts(o, ev, centIdx, ocfg, leadIsoTight, leadTruth);

        if (leadIsoTight)
        {
          fillRecoJetReplay(o, ev, centIdx, ocfg, *leadIsoTight, evJets, false);
        }
        if (leadIsoNonTight)
        {
          fillRecoJetReplay(o, ev, centIdx, ocfg, *leadIsoNonTight, evJets, true);
        }
        if (leadTruth)
        {
          fillTruthJetReplay(o, ev, centIdx, ocfg, *leadTruth, evJets);
        }
        fillSimUnfoldResponseReplay(o, ev, centIdx, ocfg, leadIsoTight, leadTruth, evJets);
      }
      o.activeView = nullptr;
    }
  }

  for (Output& o : outputs)
  {
    o.file = TFile::Open(o.outRoot.c_str(), "RECREATE");
    if (!o.file || o.file->IsZombie())
    {
      std::cerr << "[RecoilJetsPoolReplay] Could not open output " << o.outRoot << "\n";
      continue;
    }

    for (auto& scopePair : o.h)
    {
      TDirectory* dir = ensurePath(*o.file, scopePair.first);
      if (dir) dir->cd();
      for (auto& objPair : scopePair.second)
      {
        if (objPair.second) objPair.second->Write(objPair.first.c_str(), TObject::kOverwrite);
      }
    }

    if (!o.views.empty() && !o.views.front().cfg.yamlText.empty())
    {
      o.file->cd();
      TObjString yaml(o.views.front().cfg.yamlText.c_str());
      yaml.Write("analysis_config_yaml", TObject::kOverwrite);
    }
    o.file->cd();
    if (writeViewCatalog)
    {
      std::string viewId;
      std::string cfgTag;
      std::string preselection;
      std::string tight;
      std::string nonTight;
      std::string materialize;
      std::string histBase;
      std::string yamlPath;
      float jetPtMin = 0.0f;
      float dphiMin = 0.0f;
      float dphiMinPiFraction = 0.0f;
      float vzCut = 0.0f;
      float coneR = 0.0f;
      float isoFixed = 0.0f;
      int slidingIso = 0;
      int legacyTopLevel = 0;
      TTree catalog("AnalysisViewCatalog", "Catalog of replay/materialized analysis views in this ROOT file");
      catalog.Branch("view_id", &viewId);
      catalog.Branch("cfg_tag", &cfgTag);
      catalog.Branch("preselection", &preselection);
      catalog.Branch("tight", &tight);
      catalog.Branch("nonTight", &nonTight);
      catalog.Branch("materialize", &materialize);
      catalog.Branch("hist_base", &histBase);
      catalog.Branch("yaml_path", &yamlPath);
      catalog.Branch("jet_pt_min", &jetPtMin, "jet_pt_min/F");
      catalog.Branch("dphi_min", &dphiMin, "dphi_min/F");
      catalog.Branch("dphi_min_pi_fraction", &dphiMinPiFraction, "dphi_min_pi_fraction/F");
      catalog.Branch("vz_cut_cm", &vzCut, "vz_cut_cm/F");
      catalog.Branch("coneR", &coneR, "coneR/F");
      catalog.Branch("isSlidingIso", &slidingIso, "isSlidingIso/I");
      catalog.Branch("fixedGeV", &isoFixed, "fixedGeV/F");
      catalog.Branch("legacy_top_level", &legacyTopLevel, "legacy_top_level/I");
      for (const View& view : o.views)
      {
        viewId = view.viewId;
        cfgTag = view.entry.cfgTag;
        preselection = view.entry.pre;
        tight = view.entry.tight;
        nonTight = view.entry.nonTight;
        materialize = view.materialize;
        histBase = view.legacyTopLevel ? std::string{} : ("views/" + view.viewId);
        yamlPath = view.entry.yamlPath;
        jetPtMin = static_cast<float>(view.cfg.minJetPt);
        dphiMin = static_cast<float>(view.cfg.minBackToBack);
        dphiMinPiFraction = static_cast<float>(view.cfg.minBackToBack / M_PI);
        vzCut = static_cast<float>(view.cfg.vzCutCm);
        coneR = static_cast<float>(view.cfg.coneR);
        slidingIso = view.cfg.slidingIso ? 1 : 0;
        isoFixed = static_cast<float>(view.cfg.isoFixed);
        legacyTopLevel = view.legacyTopLevel ? 1 : 0;
        catalog.Fill();
        if (!view.entry.yamlPath.empty() && !isInlineViewSpec(view.entry.yamlPath) && !view.cfg.yamlText.empty())
        {
          TObjString viewYaml(view.cfg.yamlText.c_str());
          const std::string yamlName = "analysis_config_yaml_" + view.viewId;
          viewYaml.Write(yamlName.c_str(), TObject::kOverwrite);
        }
      }
      catalog.Write("", TObject::kOverwrite);
    }
    o.file->Write("", TObject::kOverwrite);
    o.file->Close();
    delete o.file;
    o.file = nullptr;
  }
}

namespace rjpool
{
  struct ObjSig
  {
    std::string cls;
    int dim = 0;
    int nx = 0, ny = 0, nz = 0;
    double xlo = 0.0, xhi = 0.0;
    double ylo = 0.0, yhi = 0.0;
    double zlo = 0.0, zhi = 0.0;
  };

  ObjSig makeObjSig(TObject* obj)
  {
    ObjSig s;
    if (!obj) return s;
    s.cls = obj->ClassName();
    if (auto* h = dynamic_cast<TH1*>(obj))
    {
      s.dim = h->GetDimension();
      s.nx = h->GetXaxis() ? h->GetXaxis()->GetNbins() : 0;
      s.ny = h->GetYaxis() ? h->GetYaxis()->GetNbins() : 0;
      s.nz = h->GetZaxis() ? h->GetZaxis()->GetNbins() : 0;
      s.xlo = h->GetXaxis() ? h->GetXaxis()->GetXmin() : 0.0;
      s.xhi = h->GetXaxis() ? h->GetXaxis()->GetXmax() : 0.0;
      s.ylo = h->GetYaxis() ? h->GetYaxis()->GetXmin() : 0.0;
      s.yhi = h->GetYaxis() ? h->GetYaxis()->GetXmax() : 0.0;
      s.zlo = h->GetZaxis() ? h->GetZaxis()->GetXmin() : 0.0;
      s.zhi = h->GetZaxis() ? h->GetZaxis()->GetXmax() : 0.0;
    }
    return s;
  }

  bool sameObjSig(const ObjSig& a, const ObjSig& b)
  {
    auto close = [](double x, double y) { return std::fabs(x - y) < 1e-9; };
    return a.cls == b.cls && a.dim == b.dim &&
           a.nx == b.nx && a.ny == b.ny && a.nz == b.nz &&
           close(a.xlo, b.xlo) && close(a.xhi, b.xhi) &&
           close(a.ylo, b.ylo) && close(a.yhi, b.yhi) &&
           close(a.zlo, b.zlo) && close(a.zhi, b.zhi);
  }

  void collectObjSigs(TDirectory* dir, const std::string& prefix, std::map<std::string, ObjSig>& out)
  {
    if (!dir) return;
    TIter next(dir->GetListOfKeys());
    while (auto* keyObj = next())
    {
      auto* key = dynamic_cast<TKey*>(keyObj);
      if (!key) continue;
      std::unique_ptr<TObject> obj(key->ReadObj());
      if (!obj) continue;
      const std::string path = prefix.empty() ? key->GetName() : prefix + "/" + key->GetName();
      if (obj->InheritsFrom(TDirectory::Class()))
      {
        collectObjSigs(dynamic_cast<TDirectory*>(obj.get()), path, out);
      }
      else
      {
        out[path] = makeObjSig(obj.get());
      }
    }
  }

  double histIntegralForPath(TFile& f, const std::string& path)
  {
    auto* obj = f.Get(path.c_str());
    auto* h = dynamic_cast<TH1*>(obj);
    return h ? h->Integral() : std::numeric_limits<double>::quiet_NaN();
  }
}

void CompareRecoilJetsLegacyParity(const char* directRoot,
                                   const char* replayRoot,
                                   const char* reportPath = "legacy_output_parity_report.txt")
{
  using namespace rjpool;
  TFile direct(directRoot, "READ");
  TFile replay(replayRoot, "READ");
  std::ofstream report(reportPath ? reportPath : "legacy_output_parity_report.txt");
  if (direct.IsZombie() || replay.IsZombie())
  {
    std::cerr << "legacy_output_parity=FAIL reason=cannot_open_input\n";
    if (report) report << "legacy_output_parity=FAIL\nreason=cannot_open_input\n";
    return;
  }

  std::map<std::string, ObjSig> dcat, rcat;
  collectObjSigs(&direct, "", dcat);
  collectObjSigs(&replay, "", rcat);

  std::vector<std::string> missing;
  std::vector<std::string> extra;
  std::vector<std::string> mismatched;
  for (const auto& [path, sig] : dcat)
  {
    auto it = rcat.find(path);
    if (it == rcat.end()) missing.push_back(path);
    else if (!sameObjSig(sig, it->second)) mismatched.push_back(path);
  }
  for (const auto& [path, sig] : rcat)
  {
    if (dcat.find(path) == dcat.end()) extra.push_back(path);
  }

  const bool pass = missing.empty() && mismatched.empty();
  std::ostream& os = report ? report : std::cout;
  os << "legacy_output_parity=" << (pass ? "PASS" : "FAIL") << "\n";
  os << "direct=" << directRoot << "\n";
  os << "replay=" << replayRoot << "\n";
  os << "direct_object_count=" << dcat.size() << "\n";
  os << "replay_object_count=" << rcat.size() << "\n";
  os << "missing_count=" << missing.size() << "\n";
  os << "mismatched_class_or_binning_count=" << mismatched.size() << "\n";
  os << "extra_count=" << extra.size() << "\n";

  auto writeList = [&](const char* header, const std::vector<std::string>& xs)
  {
    os << header << "\n";
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
      if (i >= 500)
      {
        os << "  ... truncated at 500\n";
        break;
      }
      os << "  " << xs[i] << "\n";
    }
  };
  writeList("missing_objects:", missing);
  writeList("mismatched_class_or_binning:", mismatched);
  writeList("extra_replay_objects:", extra);

  static const std::vector<std::string> probes = {
    "SIM/h_Eiso_pT_15_17",
    "SIM/h_Eiso_tight_pT_15_17",
    "SIM/h_isIsolated_isTight_pT_15_17",
    "SIM/h_sigABCD_MC_pT_15_17",
    "SIM/h2_unfoldReco_pTgamma_xJ_incl_r04",
    "SIM/h2_unfoldTruth_pTgamma_xJ_incl_r04",
    "SIM/h2_unfoldResponse_pTgamma_xJ_incl_r04",
    "SIM/h_JES3_pT_xJ_alpha_r04",
    "SIM/h_JES3TruthPure_pT_xJ_alpha_r04",
    "SIM/h_jetPt_all_r04"
  };
  os << "representative_integrals:\n";
  for (const std::string& p : probes)
  {
    const double di = histIntegralForPath(direct, p);
    const double ri = histIntegralForPath(replay, p);
    if (std::isfinite(di) || std::isfinite(ri))
      os << "  " << p << " direct=" << di << " replay=" << ri << " diff=" << (ri - di) << "\n";
  }

  std::cout << "legacy_output_parity=" << (pass ? "PASS" : "FAIL")
            << " missing=" << missing.size()
            << " mismatched=" << mismatched.size()
            << " report=" << (reportPath ? reportPath : "legacy_output_parity_report.txt")
            << "\n";
}
