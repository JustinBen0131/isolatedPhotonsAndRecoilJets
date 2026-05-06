#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TObjString.h>
#include <TProfile.h>
#include <TTree.h>
#include <TVector2.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
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
    float pt = std::numeric_limits<float>::quiet_NaN();
    float eta = std::numeric_limits<float>::quiet_NaN();
    float phi = std::numeric_limits<float>::quiet_NaN();
    float energy = std::numeric_limits<float>::quiet_NaN();
    float eiso = std::numeric_limits<float>::quiet_NaN();
    float eisoR03 = std::numeric_limits<float>::quiet_NaN();
    float eisoR04 = std::numeric_limits<float>::quiet_NaN();
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
  };

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

  TH1F* abcdHist(Output& out, const std::string& trig, const std::string& base, const std::string& sfx)
  {
    return hist1(out, trig, base + sfx, base + sfx + ";count;entries", 1, 0.5, 1.5);
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
      if (leadRecoTight)
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
      if (leadRecoTight && leadTruth)
      {
        if (auto* h = hist2Var(out, trig, "h2_unfoldResponsePho_pTgamma" + cs,
                               "h2_unfoldResponsePho_pTgamma" + cs + ";p_{T}^{#gamma,truth} [GeV];p_{T}^{#gamma,reco} [GeV]",
                               cfg.unfoldTruthPhotonPtBins, cfg.unfoldRecoPhotonPtBins))
          h->Fill(leadTruth->pt, leadRecoTight->pt, ev.weight);
      }
      else if (leadRecoTight)
      {
        if (auto* h = hist1Var(out, trig, "h_unfoldRecoPhoFakes_pTgamma" + cs,
                               "h_unfoldRecoPhoFakes_pTgamma" + cs + ";p_{T}^{#gamma,reco} [GeV];Reco fakes",
                               cfg.unfoldRecoPhotonPtBins))
          h->Fill(leadRecoTight->pt, ev.weight);
      }
      else if (leadTruth)
      {
        if (auto* h = hist1Var(out, trig, "h_unfoldTruthPhoMisses_pTgamma" + cs,
                               "h_unfoldTruthPhoMisses_pTgamma" + cs + ";p_{T}^{#gamma,truth} [GeV];Truth misses",
                               cfg.unfoldTruthPhotonPtBins))
          h->Fill(leadTruth->pt, ev.weight);
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
        const double overlapR = std::max(0.30, jetRFromKey(rKey));
        for (const JetRec* j : fidJets)
        {
          if (!j || j == recoil1) continue;
          if (dR(j->eta, j->phi, leadPho.eta, leadPho.phi) <= overlapR) continue;
          if (j->pt > jet2Pt) jet2Pt = j->pt;
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
          if (auto* h = hist3VarUniform(out, trig, "h_JES3_pT_xJ_alpha_" + rKey + cs,
                                        "h_JES3_pT_xJ_alpha_" + rKey + cs + ";p_{T}^{#gamma} [GeV];x_{J};#alpha",
                                        cfg.ptBins, 60, 0.0, 3.0, 40, 0.0, 2.0))
            h->Fill(leadPho.pt, xJ, alpha, ev.weight);
          if (auto* h = hist3VarUniform(out, trig, "h_JES3_pT_jet1Pt_alpha_" + rKey + cs,
                                        "h_JES3_pT_jet1Pt_alpha_" + rKey + cs + ";p_{T}^{#gamma} [GeV];p_{T}^{jet1} [GeV];#alpha",
                                        cfg.ptBins, 120, 0.0, 60.0, 40, 0.0, 2.0))
            h->Fill(leadPho.pt, recoil1->pt, alpha, ev.weight);
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
  if (!gForcedFanoutEntries && uniqueRoots.size() > maxOpenOutputs)
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
                     {"schema", "run", "evt", "eventKey", "isAuAu", "centBin", "centPercent", "vz", "weight", "triggers"}, true) &&
    validateBranches(photons, "AnalysisPhotonPool",
                     {"schema", "eventKey", "pt", "eta", "phi", "energy", "eiso", "eiso_r03", "eiso_r04",
                      "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx", "weta35_cogx", "wphi53_cogx",
                      "et1", "et2", "et3", "et4", "e11_over_e33", "e32_over_e35",
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
      o.file = TFile::Open(e.outRoot.c_str(), "RECREATE");
      if (!o.file || o.file->IsZombie())
      {
        std::cerr << "[RecoilJetsPoolReplay] Could not open output " << e.outRoot << "\n";
        continue;
      }
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
      v.legacyTopLevel = singleView || v.viewId == "legacy" || v.viewId == "nominal_legacy";
    }
  }

  int eventRun = 0;
  int eventEvt = 0;
  long long eventKey = 0;
  int isAuAu = 0;
  int centBin = -1;
  float centPercent = -1.0f;
  float vz = 0.0f;
  float weight = 1.0f;
  std::vector<std::string>* triggers = nullptr;
  events.SetBranchAddress("run", &eventRun);
  events.SetBranchAddress("evt", &eventEvt);
  events.SetBranchAddress("eventKey", &eventKey);
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
  float weta = 0, wphi = 0, weta33 = 0, wphi33 = 0, weta35 = 0, wphi53 = 0;
  float et1 = 0, et2 = 0, et3 = 0, et4 = 0, e11e33 = 0, e32e35 = 0;
  float e11e22 = 0, e11e13 = 0, e11e15 = 0, e11e17 = 0, e11e31 = 0, e11e51 = 0, e11e71 = 0;
  float e22e33 = 0, e22e35 = 0, e22e37 = 0, e22e53 = 0, w32 = 0, w52 = 0, w72 = 0, meanTime = 0;
  float npb = 0, tightBDT = 0, auauNpb = 0, auauTightBDT = 0;
  float mbdTime = 0, clusterMbdDeltaT = 0;
  int npbHasAwayJet = 0, npbLabel = -1, isNPB = -1;
  int truthSignal = 0, truthTrackId = -1;
  float truthPt = 0, truthEta = 0, truthPhi = 0;
  photons.SetBranchAddress("eventKey", &eventKey);
  photons.SetBranchAddress("pt", &pt);
  photons.SetBranchAddress("eta", &eta);
  photons.SetBranchAddress("phi", &phi);
  photons.SetBranchAddress("energy", &energy);
  photons.SetBranchAddress("eiso", &eiso);
  photons.SetBranchAddress("eiso_r03", &eisoR03);
  photons.SetBranchAddress("eiso_r04", &eisoR04);
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

  const Long64_t nPho = photons.GetEntries();
  for (Long64_t i = 0; i < nPho; ++i)
  {
    photons.GetEntry(i);
    const long long key = scopedEventKey(photons, eventKey);
    const auto evIt = eventMap.find(key);
    if (evIt == eventMap.end()) continue;
    PhotonRec pr;
    pr.eventKey = key;
    pr.pt = pt;
    pr.eta = eta;
    pr.phi = phi;
    pr.energy = energy;
    pr.eiso = eiso;
    pr.eisoR03 = eisoR03;
    pr.eisoR04 = eisoR04;
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
    photonsByEvent[key].push_back(pr);
  }

  std::unordered_map<long long, std::vector<JetRec>> jetsByEvent;
  long long jetEventKey = 0;
  int jetIsTruth = 0;
  float jetPt = 0, jetRawPt = 0, jetAreaSubPt = 0, jetEta = 0, jetPhi = 0, jetArea = 0;
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

        const PhotonRec* leadIsoTight = nullptr;
        const PhotonRec* leadIsoNonTight = nullptr;

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
          if (!iso && !nonIso) continue;

          if (!passPre(view.entry, ocfg, ph.pt, ph.weta, ph.et1, ph.e11e33, ph.e32e35, ph.npb, ph.auauNpb)) continue;
          const int tag = tightTag(view.entry, ocfg, ph.pt, ph.weta, ph.wphi, ph.et1, ph.e11e33, ph.e32e35, ph.tightBDT, ph.auauTightBDT);
          if (tag == 3) continue;

          const char* base = nullptr;
          int regionBin = 0;
          if (iso && tag == 1) { base = "h_isIsolated_isTight"; regionBin = 1; }
          else if (nonIso && tag == 1) { base = "h_notIsolated_isTight"; regionBin = 2; }
          else if (iso && tag == 2) { base = "h_isIsolated_notTight"; regionBin = 3; }
          else if (nonIso && tag == 2) { base = "h_notIsolated_notTight"; regionBin = 4; }
          if (!base) continue;

          const std::string sfx = suffix(ptIdx, centIdx, ev.isAuAu, ocfg);
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
        }

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
      }
      o.activeView = nullptr;
    }
  }

  for (Output& o : outputs)
  {
    if (!o.file) continue;
    if (!o.views.empty() && !o.views.front().cfg.yamlText.empty())
    {
      o.file->cd();
      TObjString yaml(o.views.front().cfg.yamlText.c_str());
      yaml.Write("analysis_config_yaml", TObject::kOverwrite);
    }
    o.file->cd();
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
    o.file->Write("", TObject::kOverwrite);
    o.file->Close();
    delete o.file;
    o.file = nullptr;
  }
}
