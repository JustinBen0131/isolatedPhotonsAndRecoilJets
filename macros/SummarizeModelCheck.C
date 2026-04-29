#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct HistTotals {
  double entries = 0.0;
  double integral = 0.0;
  double weightedMean = 0.0;
  double underflow = 0.0;
  double overflow = 0.0;
  double bins[8] = {0.0};
  double minFilled = std::numeric_limits<double>::infinity();
  double maxFilled = -std::numeric_limits<double>::infinity();
  double sentinelLike = 0.0;
  double tightWindow = 0.0;
  double npbPassWindow = 0.0;
  double aboveZero = 0.0;
  int histCount = 0;
};

struct SummaryRow {
  std::string cfg;
  std::string label;
  double a = 0.0;
  double b = 0.0;
  double c = 0.0;
  double d = 0.0;
  double cOverA = -1.0;
  double sigA = 0.0;
  double sigB = 0.0;
  double sigC = 0.0;
  double sigD = 0.0;
  double sigCOverA = -1.0;
  double sigCOverAC = -1.0;
  double sigCOverRecoC = -1.0;
  double tightFinite = 0.0;
  double tightMissing = 0.0;
  double tightMissFrac = -1.0;
  HistTotals tightAllScore;
  HistTotals tightPreScore;
  double npbFinite = 0.0;
  double npbMissing = 0.0;
  double npbMissFrac = -1.0;
  HistTotals npbAllScore;
};

std::string trim(const std::string &s)
{
  const char *ws = " \t\r\n";
  const size_t b = s.find_first_not_of(ws);
  if (b == std::string::npos) return "";
  const size_t e = s.find_last_not_of(ws);
  return s.substr(b, e - b + 1);
}

std::string shellQuote(const std::string &s)
{
  std::string out = "'";
  for (char c : s) {
    if (c == '\'') out += "'\\''";
    else out += c;
  }
  out += "'";
  return out;
}

bool startsWith(const std::string &s, const std::string &prefix)
{
  return s.rfind(prefix, 0) == 0;
}

bool endsWith(const std::string &s, const std::string &suffix)
{
  return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::string> readInputFiles(const std::string &input)
{
  std::vector<std::string> files;

  if (endsWith(input, ".root")) {
    files.push_back(input);
    return files;
  }

  if (endsWith(input, ".txt") || endsWith(input, ".list")) {
    std::ifstream in(input);
    std::string line;
    while (std::getline(in, line)) {
      line = trim(line);
      if (line.empty() || line[0] == '#') continue;
      files.push_back(line);
    }
    return files;
  }

  const std::string findCmd = "find " + shellQuote(input) + " -type f -name '*.root' | sort";
  TString found = gSystem->GetFromPipe(findCmd.c_str());
  std::unique_ptr<TObjArray> lines(found.Tokenize("\n"));
  for (int i = 0; i < lines->GetEntriesFast(); ++i) {
    const auto *os = static_cast<TObjString *>(lines->At(i));
    if (!os) continue;
    const std::string path = trim(os->GetString().Data());
    if (!path.empty()) files.push_back(path);
  }
  return files;
}

void accumulateMatching(TDirectory *dir, const std::string &prefix, HistTotals &totals)
{
  if (!dir) return;

  TIter next(dir->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next())) {
    TObject *obj = key->ReadObj();
    if (!obj) continue;

    if (obj->InheritsFrom(TDirectory::Class())) {
      accumulateMatching(static_cast<TDirectory *>(obj), prefix, totals);
    } else if (obj->InheritsFrom(TH1::Class())) {
      TH1 *h = static_cast<TH1 *>(obj);
      const std::string name = h->GetName();
      if (!startsWith(name, prefix)) continue;

      ++totals.histCount;
      const double entries = h->GetEntries();
      totals.entries += entries;
      totals.integral += h->Integral(0, h->GetNbinsX() + 1);
      totals.weightedMean += h->GetMean() * entries;
      totals.underflow += h->GetBinContent(0);
      totals.overflow += h->GetBinContent(h->GetNbinsX() + 1);

      const int nbAbcd = std::min(h->GetNbinsX(), 7);
      for (int ib = 1; ib <= nbAbcd; ++ib) {
        totals.bins[ib] += h->GetBinContent(ib);
      }

      for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        if (h->GetBinContent(ib) <= 0.0) continue;
        const double x = h->GetXaxis()->GetBinCenter(ib);
        const double y = h->GetBinContent(ib);
        totals.minFilled = std::min(totals.minFilled, x);
        totals.maxFilled = std::max(totals.maxFilled, x);
        if (x < -0.95) totals.sentinelLike += y;
        if (x > 0.8 && x < 1.0) totals.tightWindow += y;
        if (x > 0.5) totals.npbPassWindow += y;
        if (x > 0.0) totals.aboveZero += y;
      }
    }
  }
}

std::string configTagFromPath(const std::string &path)
{
  const std::string needle = "/model_check/";
  const size_t pos = path.find(needle);
  if (pos != std::string::npos) {
    const size_t start = pos + needle.size();
    const size_t slash = path.find('/', start);
    if (slash != std::string::npos) return path.substr(start, slash - start);
  }

  const size_t slash = path.find_last_of('/');
  const std::string file = (slash == std::string::npos) ? path : path.substr(slash + 1);
  const std::string prefix = "RecoilJets_isSim_sim_run28_photonjet20_";
  const std::string suffix = "_CHECKMODELS_firstfile_grp001.root";
  if (startsWith(file, prefix) && endsWith(file, suffix)) {
    return file.substr(prefix.size(), file.size() - prefix.size() - suffix.size());
  }
  return file;
}

std::string photonIdLabel(const std::string &cfg)
{
  if (cfg.find("preselectionReference_tightReference_nonTightReference") != std::string::npos) {
    return "reference";
  }
  if (cfg.find("preselectionVariantA_tightVariantA_nonTightVariantA") != std::string::npos ||
      cfg.find("preselectionVariantA_tightVariantA_nonTightReference") != std::string::npos) {
    return "variantA pre + tight/non-tight BDT";
  }
  if (cfg.find("preselectionVariantB_tightVariantA_nonTightVariantA") != std::string::npos ||
      cfg.find("preselectionVariantB_tightVariantA_nonTightReference") != std::string::npos) {
    return "variantB pre + tight/non-tight BDT";
  }
  return "unknown";
}

std::string passWarn(bool pass)
{
  return pass ? "PASS" : "WARN";
}

double meanOf(const HistTotals &h)
{
  return h.entries > 0.0 ? h.weightedMean / h.entries : -999.0;
}

double minOf(const HistTotals &h)
{
  return std::isfinite(h.minFilled) ? h.minFilled : -999.0;
}

double maxOf(const HistTotals &h)
{
  return std::isfinite(h.maxFilled) ? h.maxFilled : -999.0;
}

double fraction(double num, double den)
{
  return den > 0.0 ? num / den : -1.0;
}

void printScoreBlock(const std::string &name, const HistTotals &h)
{
  std::cout << "    " << name
            << ": entries=" << h.entries
            << " hists=" << h.histCount
            << " mean=" << meanOf(h)
            << " minBin=" << minOf(h)
            << " maxBin=" << maxOf(h)
            << " sentinelLike=" << h.sentinelLike
            << " tightWindow(0.8,1.0)=" << h.tightWindow
            << " npbPass(>0.5)=" << h.npbPassWindow
            << " positive=" << h.aboveZero
            << " underflow=" << h.underflow
            << " overflow=" << h.overflow
            << "\n";
}

void printQuickScoreVerdict(const SummaryRow &row)
{
  const bool usesTightBdt = row.cfg.find("tightVariantA") != std::string::npos;
  const bool usesNpb = row.cfg.find("preselectionVariantA") != std::string::npos;

  if (usesTightBdt) {
    const double sentinelFrac = fraction(row.tightAllScore.sentinelLike, row.tightAllScore.entries);
    const bool tightScoreHealthy =
      row.tightAllScore.entries > 0.0 &&
      row.tightMissing == 0.0 &&
      row.tightAllScore.tightWindow > 0.0 &&
      sentinelFrac < 0.5;
    std::cout << "    QUICK tight-BDT check: "
              << "entries=" << row.tightAllScore.entries
              << " sentinelFrac=" << sentinelFrac
              << " tightWindowAll=" << row.tightAllScore.tightWindow
              << " tightWindowPre=" << row.tightPreScore.tightWindow
              << " maxBin=" << maxOf(row.tightAllScore)
              << " " << passWarn(tightScoreHealthy)
              << "\n";
  } else {
    std::cout << "    QUICK tight-BDT check: not enabled for this reference row PASS\n";
  }

  if (usesNpb) {
    const bool npbHealthy =
      row.npbAllScore.entries > 0.0 &&
      row.npbMissing == 0.0 &&
      row.npbAllScore.npbPassWindow > 0.0;
    std::cout << "    QUICK NPB check: "
              << "entries=" << row.npbAllScore.entries
              << " npbPass(>0.5)=" << row.npbAllScore.npbPassWindow
              << " positive=" << row.npbAllScore.aboveZero
              << " maxBin=" << maxOf(row.npbAllScore)
              << " " << passWarn(npbHealthy)
              << "\n";
  } else {
    std::cout << "    QUICK NPB check: not enabled for this row PASS\n";
  }
}

void printComparison(const SummaryRow &ref, const SummaryRow *row, const char *name)
{
  if (!row) {
    std::cout << "  " << name << " missing WARN\n";
    return;
  }

  const double dC = ref.c - row->c;
  const double dSigC = ref.sigC - row->sigC;
  const double fracSigC = ref.sigC > 0.0 ? dSigC / ref.sigC : -1.0;
  std::cout << "  " << name
            << ": recoC=" << row->c
            << " recoC/A=" << row->cOverA
            << " sigA=" << row->sigA
            << " sigC=" << row->sigC
            << " sigC/A=" << row->sigCOverA
            << " sigC/(A+C)=" << row->sigCOverAC
            << " sigC/recoC=" << row->sigCOverRecoC
            << " deltaRecoC(ref-current)=" << dC
            << " deltaSigC(ref-current)=" << dSigC
            << " fracSigCReduction=" << fracSigC
            << " " << passWarn(row->sigC <= ref.sigC || ref.sigC <= 0.0)
            << "\n";
}

}  // namespace

void SummarizeModelCheck(const char *inputPath = "local_sim_outputs/model_check", const char *logPath = "")
{
  const std::string input = inputPath ? inputPath : "";
  const std::string log = logPath ? logPath : "";
  if (input.empty()) {
    std::cerr << "[model-check] Empty input path" << std::endl;
    return;
  }

  const std::vector<std::string> files = readInputFiles(input);

  std::cout << "\n================ MODEL CHECK SUMMARY ================\n";
  std::cout << "Input: " << input << "\n";
  if (!log.empty()) {
    std::cout << "Summary log for sftp: " << log << "\n";
  }
  std::cout << "Files: " << files.size() << "\n\n";

  if (files.empty()) {
    std::cout << "[FAIL] No ROOT files found in input path or manifest.\n";
    std::cout << "=====================================================\n";
    return;
  }

  int nFilesOpened = 0;
  int nFilesWithABCD = 0;
  int nTightWarnings = 0;
  int nNpbWarnings = 0;
  std::vector<SummaryRow> rows;

  std::cout << std::fixed << std::setprecision(4);
  std::cout << "config"
            << " | reco A B C D"
            << " | reco C/A"
            << " | sig A B C D"
            << " | sig C/A"
            << " | sig C/(A+C)"
            << " | sig C/reco C"
            << " | tight finite/missing/missFrac"
            << " | npb finite/missing/missFrac"
            << "\n";
  std::cout << "--------------------------------------------------------------------------------\n";

  for (const std::string &path : files) {
    TFile f(path.c_str(), "READ");
    if (f.IsZombie()) {
      std::cout << "[WARN] Could not open " << path << "\n";
      continue;
    }
    ++nFilesOpened;

    HistTotals hA, hB, hC, hD;
    HistTotals tightFinite, tightMissing, tightAllScore, tightPreScore;
    HistTotals npbFinite, npbMissing, npbAllScore;
    HistTotals sigAbcd;

    accumulateMatching(&f, "h_pTgamma_ABCD_A", hA);
    accumulateMatching(&f, "h_pTgamma_ABCD_B", hB);
    accumulateMatching(&f, "h_pTgamma_ABCD_C", hC);
    accumulateMatching(&f, "h_pTgamma_ABCD_D", hD);
    accumulateMatching(&f, "h_tightBDTScore_finite", tightFinite);
    accumulateMatching(&f, "h_tightBDTScore_missing", tightMissing);
    accumulateMatching(&f, "h_tightBDTScore_allCandidates", tightAllScore);
    accumulateMatching(&f, "h_tightBDTScore_preselected", tightPreScore);
    accumulateMatching(&f, "h_npbScore_finite", npbFinite);
    accumulateMatching(&f, "h_npbScore_missing", npbMissing);
    accumulateMatching(&f, "h_npbScore_allCandidates", npbAllScore);
    accumulateMatching(&f, "h_sigABCD_MC", sigAbcd);

    SummaryRow row;
    row.cfg = configTagFromPath(path);
    row.label = photonIdLabel(row.cfg);
    row.a = hA.integral;
    row.b = hB.integral;
    row.c = hC.integral;
    row.d = hD.integral;
    row.cOverA = row.a > 0.0 ? row.c / row.a : -1.0;
    row.sigA = sigAbcd.bins[1];
    row.sigB = sigAbcd.bins[2];
    row.sigC = sigAbcd.bins[3];
    row.sigD = sigAbcd.bins[4];
    row.sigCOverA = row.sigA > 0.0 ? row.sigC / row.sigA : -1.0;
    row.sigCOverAC = (row.sigA + row.sigC) > 0.0 ? row.sigC / (row.sigA + row.sigC) : -1.0;
    row.sigCOverRecoC = row.c > 0.0 ? row.sigC / row.c : -1.0;
    row.tightFinite = tightFinite.entries;
    row.tightMissing = tightMissing.entries;
    row.tightMissFrac = (row.tightFinite + row.tightMissing) > 0.0 ? row.tightMissing / (row.tightFinite + row.tightMissing) : -1.0;
    row.tightAllScore = tightAllScore;
    row.tightPreScore = tightPreScore;
    row.npbFinite = npbFinite.entries;
    row.npbMissing = npbMissing.entries;
    row.npbMissFrac = (row.npbFinite + row.npbMissing) > 0.0 ? row.npbMissing / (row.npbFinite + row.npbMissing) : -1.0;
    row.npbAllScore = npbAllScore;

    if ((row.a + row.b + row.c + row.d) > 0.0) ++nFilesWithABCD;

    const bool usesTightBdt = row.cfg.find("tightVariantA") != std::string::npos;
    const bool tightOk = !usesTightBdt ||
      (row.tightFinite > 0.0 &&
       row.tightMissing == 0.0 &&
       row.tightAllScore.entries > 0.0 &&
       row.tightAllScore.tightWindow > 0.0 &&
       fraction(row.tightAllScore.sentinelLike, row.tightAllScore.entries) < 0.5);
    const bool usesNpb = row.cfg.find("preselectionVariantA") != std::string::npos;
    const bool npbOk = !usesNpb ||
      (row.npbFinite > 0.0 &&
       row.npbMissing == 0.0 &&
       row.npbAllScore.entries > 0.0 &&
       row.npbAllScore.npbPassWindow > 0.0);
    if (!tightOk) ++nTightWarnings;
    if (!npbOk) ++nNpbWarnings;

    rows.push_back(row);

    std::cout << row.cfg
              << " | " << row.a << " " << row.b << " " << row.c << " " << row.d
              << " | " << row.cOverA
              << " | " << row.sigA << " " << row.sigB << " " << row.sigC << " " << row.sigD
              << " | " << row.sigCOverA
              << " | " << row.sigCOverAC
              << " | " << row.sigCOverRecoC
              << " | " << row.tightFinite << "/" << row.tightMissing << "/" << row.tightMissFrac
              << " " << passWarn(tightOk)
              << " | " << row.npbFinite << "/" << row.npbMissing << "/" << row.npbMissFrac
              << " " << passWarn(npbOk)
              << "\n";
    printScoreBlock("tight score all candidates", row.tightAllScore);
    printScoreBlock("tight score after preselection", row.tightPreScore);
    printScoreBlock("NPB score all candidates", row.npbAllScore);
    printQuickScoreVerdict(row);
  }

  std::cout << "\nChecklist:\n";
  std::cout << "  opened files          : " << nFilesOpened << " / " << files.size()
            << " " << passWarn(nFilesOpened == static_cast<int>(files.size())) << "\n";
  std::cout << "  files with ABCD counts: " << nFilesWithABCD << " / " << nFilesOpened
            << " " << passWarn(nFilesWithABCD == nFilesOpened) << "\n";
  std::cout << "  tight BDT diagnostics : warnings=" << nTightWarnings
            << " " << passWarn(nTightWarnings == 0) << "\n";
  std::cout << "  NPB diagnostics       : warnings=" << nNpbWarnings
            << " " << passWarn(nNpbWarnings == 0) << "\n";

  const SummaryRow *reference = nullptr;
  const SummaryRow *variantA = nullptr;
  const SummaryRow *variantB = nullptr;
  for (const auto &row : rows) {
    if (row.label == "reference") reference = &row;
    else if (row.label == "variantA pre + tight/non-tight BDT") variantA = &row;
    else if (row.label == "variantB pre + tight/non-tight BDT") variantB = &row;
  }

  std::cout << "\nPhoton-ID comparison against reference:\n";
  if (!reference) {
    std::cout << "  reference row missing: cannot compute improvement checks WARN\n";
  } else {
    std::cout << "  reference"
              << ": recoC=" << reference->c
              << " recoC/A=" << reference->cOverA
              << " sigA=" << reference->sigA
              << " sigC=" << reference->sigC
              << " sigC/A=" << reference->sigCOverA
              << " sigC/(A+C)=" << reference->sigCOverAC
              << " sigC/recoC=" << reference->sigCOverRecoC
              << "\n";

    if (variantB) {
      printComparison(*reference, variantB, "variantB pre + tight/non-tight BDT");
    }
    printComparison(*reference, variantA, "variantA pre + tight/non-tight BDT");

    if (variantA && variantB) {
      const bool orderingOk = (reference->sigC <= 0.0) || (variantA->sigC <= variantB->sigC);
      std::cout << "  expected ordering variantA sigC <= variantB sigC: "
                << passWarn(orderingOk)
                << " (variantA=" << variantA->sigC
                << ", variantB=" << variantB->sigC << ")\n";
    }
  }

  std::cout << "\nInterpretation notes:\n";
  std::cout << "  h_sigABCD_MC bins are matched truth-signal reco classifications: 1=A, 2=B, 3=C, 4=D.\n";
  std::cout << "  The plotted leakage quantity is sigC/sigA. If sigA is zero, the BDT is classifying matched truth signal as non-tight.\n";
  std::cout << "  Score means/ranges above come from h_*Score_allCandidates and h_*Score_preselected, not the finite-count histograms.\n";
  std::cout << "  For the 1000-event checkModels pass, the critical checks are nonzero tightWindow(0.8,1.0), low sentinelFrac, and NPB npbPass(>0.5).\n";
  if (!log.empty()) {
    std::cout << "\nSFTP this summary log:\n  " << log << "\n";
  }
  std::cout << "=====================================================\n";
}
