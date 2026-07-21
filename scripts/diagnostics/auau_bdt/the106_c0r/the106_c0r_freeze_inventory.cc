#include "the106_c0r_diagnostic_writer.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct Row
{
  std::string directory;
  std::string object_name;
  std::string title;
  std::string variable;
  std::string trigger;
  std::string tag;
  int photon_pt_slice;
  int centrality_slice;
  int nbins;
  std::string xmin_bits;
  std::string xmax_bits;
  std::string x_axis_title;
  std::string y_axis_title;
  bool sumw2;
};

std::string encoded(const std::string& value)
{
  return the106::c0r::percentEncode(value);
}

int sliceIndex(const std::map<std::string, int>& slices,
               const std::string& low,
               const std::string& high)
{
  const std::map<std::string, int>::const_iterator found =
      slices.find(low + "_" + high);
  if (found == slices.end()) throw std::runtime_error("unsupported slice");
  return found->second;
}

bool rowLess(const Row& left, const Row& right)
{
  if (left.directory != right.directory) return left.directory < right.directory;
  return left.object_name < right.object_name;
}
}  // namespace

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "usage: the106_c0r_freeze_inventory armA1.root inventory.tsv\n";
    return 2;
  }

  try
  {
    const std::set<std::string> triggers = {
        "MBD_NS_geq_2_vtx_lt_150",
        "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
        "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"};
    const std::set<std::string> variables = {
        "weta", "wphi", "weta33", "wphi33", "weta35",
        "wphi53", "et1", "e11e33", "e32e35"};
    const std::map<std::string, int> pt_slices = {
        {"15_17", 0}, {"17_19", 1}, {"19_21", 2},
        {"21_23", 3}, {"23_26", 4}, {"26_35", 5}};
    const std::map<std::string, int> centrality_slices = {
        {"0_20", 0}, {"20_50", 1}, {"50_80", 2}};
    const std::regex name_pattern(
        "^h_ss_(weta|wphi|weta33|wphi33|weta35|wphi53|et1|e11e33|e32e35)_"
        "(pre|tight|nonTight)_pT_([0-9]+)_([0-9]+)_cent_([0-9]+)_([0-9]+)$");

    TFile input(argv[1], "READ");
    if (input.IsZombie()) throw std::runtime_error("cannot open direct ROOT input");

    std::vector<Row> rows;
    std::set<std::string> seen_paths;
    std::set<std::string> seen_variables;
    TIter top_keys(input.GetListOfKeys());
    while (TKey* top_key = static_cast<TKey*>(top_keys()))
    {
      TDirectory* directory = dynamic_cast<TDirectory*>(top_key->ReadObj());
      if (!directory) continue;
      const std::string directory_name = directory->GetName();
      TIter object_keys(directory->GetListOfKeys());
      while (TKey* object_key = static_cast<TKey*>(object_keys()))
      {
        const std::string object_name = object_key->GetName();
        std::smatch match;
        if (!std::regex_match(object_name, match, name_pattern)) continue;
        if (triggers.find(directory_name) == triggers.end())
          throw std::runtime_error("target raw-QA object in unsupported trigger directory");

        TH1F* histogram = dynamic_cast<TH1F*>(object_key->ReadObj());
        if (!histogram || histogram->GetDimension() != 1)
          throw std::runtime_error("target raw-QA object is not one-dimensional TH1F");
        const std::string variable = match[1].str();
        const std::string tag = match[2].str();
        if (variables.find(variable) == variables.end())
          throw std::runtime_error("unsupported raw-QA variable");
        const std::string expected_title = "h_ss_" + variable + "_" + tag;
        if (histogram->GetTitle() != expected_title ||
            std::string(histogram->GetXaxis()->GetTitle()) != variable ||
            std::string(histogram->GetYaxis()->GetTitle()) != "Entries" ||
            histogram->GetNbinsX() != 120 ||
            the106::c0r::binary64Hex(histogram->GetXaxis()->GetBinLowEdge(1)) !=
                the106::c0r::binary64Hex(0.0) ||
            the106::c0r::binary64Hex(histogram->GetXaxis()->GetBinUpEdge(120)) !=
                the106::c0r::binary64Hex(1.2))
          throw std::runtime_error("direct raw-QA object violates authoritative booking contract");

        Row row;
        row.directory = directory_name;
        row.object_name = object_name;
        row.title = histogram->GetTitle();
        row.variable = variable;
        row.trigger = directory_name;
        row.tag = tag;
        row.photon_pt_slice = sliceIndex(pt_slices, match[3].str(), match[4].str());
        row.centrality_slice =
            sliceIndex(centrality_slices, match[5].str(), match[6].str());
        row.nbins = histogram->GetNbinsX();
        row.xmin_bits =
            the106::c0r::binary64Hex(histogram->GetXaxis()->GetBinLowEdge(1));
        row.xmax_bits =
            the106::c0r::binary64Hex(histogram->GetXaxis()->GetBinUpEdge(row.nbins));
        row.x_axis_title = histogram->GetXaxis()->GetTitle();
        row.y_axis_title = histogram->GetYaxis()->GetTitle();
        row.sumw2 = histogram->GetSumw2N() != 0;
        const std::string path = row.directory + "/" + row.object_name;
        if (!seen_paths.insert(path).second)
          throw std::runtime_error("duplicate direct raw-QA object path");
        seen_variables.insert(variable);
        rows.push_back(row);
      }
    }
    if (rows.empty() || seen_variables != variables)
      throw std::runtime_error("direct bounded output lacks the complete nine-variable raw-QA family");
    std::sort(rows.begin(), rows.end(), rowLess);

    std::ofstream output(argv[2], std::ios::out | std::ios::trunc);
    if (!output) throw std::runtime_error("cannot create frozen inventory");
    output << "directory\tobject_name\tobject_class\ttitle\tvariable\ttrigger\ttag\t"
              "view_suffix\tphoton_pt_slice\tcentrality_slice\tnbins\txmin_bits\t"
              "xmax_bits\tx_axis_title\ty_axis_title\tsumw2_required\trequired\t"
              "content_comparator\n";
    for (std::vector<Row>::const_iterator row = rows.begin(); row != rows.end(); ++row)
    {
      output << encoded(row->directory) << '\t'
             << encoded(row->object_name) << "\tTH1F\t"
             << encoded(row->title) << '\t'
             << encoded(row->variable) << '\t'
             << encoded(row->trigger) << '\t'
             << encoded(row->tag) << "\tcanonical\t"
             << row->photon_pt_slice << '\t' << row->centrality_slice << '\t'
             << row->nbins << '\t' << row->xmin_bits << '\t' << row->xmax_bits << '\t'
             << encoded(row->x_axis_title) << '\t' << encoded(row->y_axis_title) << '\t'
             << (row->sumw2 ? "1" : "0")
             << "\t1\tBITWISE_SINGLE_PROCESS_UNWEIGHTED\n";
    }
    output.close();
    if (!output) throw std::runtime_error("failed while sealing frozen inventory");
    std::cout << "FROZEN_RAWQA_OBJECTS=" << rows.size() << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << "inventory freeze failed: " << error.what() << '\n';
    return 1;
  }
}
