#ifndef RJ_SHOWER_FACTORIAL_V1_H
#define RJ_SHOWER_FACTORIAL_V1_H

#include "RJReplayFoundationV1.h"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace RJShowerFactorialV1
{
using namespace RJReplayFoundationV1;

enum class EnergySource : std::int32_t { CALIBRATED_TOWERINFO=0, RAWCLUSTER_MAP_VALUE=1 };
enum class Membership : std::int32_t { FULL_GRID=0, RAWCLUSTER_OWNED=1 };

struct Definition
{
  const char* name;
  EnergySource energy_source;
  Membership rectangular_membership;
  Membership moment_membership;
  double floor_gev;
};

inline const std::array<Definition,7>& definitions()
{
  static const std::array<Definition,7> values{{
      {"H70",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::RAWCLUSTER_OWNED,0.070},
      {"H0",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::RAWCLUSTER_OWNED,0.0},
      {"G70",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::FULL_GRID,0.070},
      {"G0",EnergySource::CALIBRATED_TOWERINFO,Membership::FULL_GRID,Membership::FULL_GRID,0.0},
      {"O70",EnergySource::CALIBRATED_TOWERINFO,Membership::RAWCLUSTER_OWNED,Membership::RAWCLUSTER_OWNED,0.070},
      {"O0",EnergySource::CALIBRATED_TOWERINFO,Membership::RAWCLUSTER_OWNED,Membership::RAWCLUSTER_OWNED,0.0},
      {"R70",EnergySource::RAWCLUSTER_MAP_VALUE,Membership::RAWCLUSTER_OWNED,Membership::RAWCLUSTER_OWNED,0.070}}};
  return values;
}

inline const Definition& definition(const std::string& name)
{
  for (const auto& value : definitions()) if (name == value.name) return value;
  throw std::invalid_argument("unknown shower definition: " + name);
}

inline std::string semanticText(const Definition& value)
{
  return std::string("RJ_SHOWER_DEFINITION_FACTORIAL_V1|") + value.name +
      "|energy=" + std::to_string(static_cast<int>(value.energy_source)) +
      "|sums=" + std::to_string(static_cast<int>(value.rectangular_membership)) +
      "|moments=" + std::to_string(static_cast<int>(value.moment_membership)) +
      "|floor_gev=" + (value.floor_gev > 0.0 ? "0.070000" : "0.000000") +
      "|grid=7x7|tower_quality=TowerInfo_get_isGood_only|center_excluded_from_cogx_numerator=1";
}

inline std::string semanticSha256(const Definition& value)
{
  return sha256Hex(semanticText(value));
}

inline bool member(const ShowerCellRow& cell, Membership membership)
{
  return membership == Membership::FULL_GRID || cell.rawcluster_owned != 0;
}

inline bool selectedEnergy(const ShowerCellRow& cell,
                           const Definition& value,
                           double& energy)
{
  if (value.energy_source == EnergySource::CALIBRATED_TOWERINFO)
  {
    energy = cell.calibrated_energy;
    return cell.is_good != 0 && std::isfinite(energy) && energy > value.floor_gev;
  }
  energy = cell.rawcluster_map_value;
  return cell.rawcluster_owned != 0 && cell.rawcluster_value_present != 0 &&
      std::isfinite(energy) && energy > value.floor_gev;
}

inline ShowerFeatureViewRow buildView(const Identity128& candidateId,
                                      const Definition& value,
                                      const std::vector<ShowerCellRow>& cells,
                                      double rawCenterEta,
                                      double rawCenterPhi,
                                      const std::array<double,4>& nativeEt)
{
  if (candidateId.isNull()) throw std::invalid_argument("null shower-view candidate identity");
  if (!std::isfinite(rawCenterEta) || !std::isfinite(rawCenterPhi))
    throw std::invalid_argument("nonfinite shower center");

  const int centerEta=static_cast<int>(std::floor(rawCenterEta));
  int centerPhi=static_cast<int>(std::floor(rawCenterPhi));
  while (centerPhi<0) centerPhi+=256;
  while (centerPhi>=256) centerPhi-=256;
  const double cogEtaLocal=3.0+(rawCenterEta-std::floor(rawCenterEta)-0.5);
  const double cogPhiLocal=3.0+(rawCenterPhi-std::floor(rawCenterPhi)-0.5);
  auto deltaPhiIndex=[](int towerPhi,int referencePhi)
  {
    int delta=towerPhi-referencePhi;
    while(delta<-128)delta+=256;
    while(delta>127)delta-=256;
    return delta;
  };

  ShowerFeatureViewRow row;
  row.candidate_id=candidateId;
  row.definition_name=value.name;
  row.definition_id=makeIdentity(std::string("shower-definition|")+value.name+"|"+semanticText(value));
  row.semantic_sha256=semanticSha256(value);
  row.floor_gev=value.floor_gev;
  row.cog_eta=cogEtaLocal;
  row.cog_phi=cogPhiLocal;
  row.raw_center_eta=rawCenterEta;
  row.raw_center_phi=rawCenterPhi;
  row.center_eta_index=centerEta;
  row.center_phi_index=centerPhi;
  row.energy_source=static_cast<int>(value.energy_source);
  row.rectangular_membership=static_cast<int>(value.rectangular_membership);
  row.moment_membership=static_cast<int>(value.moment_membership);
  row.native_et1=nativeEt[0]; row.native_et2=nativeEt[1];
  row.native_et3=nativeEt[2]; row.native_et4=nativeEt[3];

  const int signPhi=cogPhiLocal>3.0?1:-1;
  for (const auto& cell : cells)
  {
    if (cell.candidate_id != candidateId) continue;
    const int i=cell.tower_eta_index-centerEta+3;
    const int j=deltaPhiIndex(cell.tower_phi_index,centerPhi)+3;
    if (i<0||i>6||j<0||j>6) continue;
    row.good_cell_count += cell.is_good != 0;
    row.owned_cell_count += cell.rawcluster_owned != 0;
    row.exact_zero_count += cell.is_good != 0 && cell.is_zero != 0;
    row.negative_count += cell.is_good != 0 && cell.is_negative != 0;
    row.nonfinite_count += cell.is_nonfinite != 0;

    double energy=std::numeric_limits<double>::quiet_NaN();
    if (!selectedEnergy(cell,value,energy)) continue;
    const bool sumMember=member(cell,value.rectangular_membership);
    const bool momentMember=member(cell,value.moment_membership);
    const int di=std::abs(i-3),dj=std::abs(j-3);
    if (sumMember)
    {
      ++row.active_sum_cell_count;
      if (i==3&&j==3) row.e11+=energy;
      if (di<=1&&dj<=1) row.e33+=energy;
      if (di<=1&&(j==3||j==3+signPhi)) row.e32+=energy;
      if (di<=1&&dj<=2) row.e35+=energy;
    }
    if (momentMember)
    {
      ++row.active_moment_cell_count;
      const double deta=static_cast<double>(i)-cogEtaLocal;
      const double dphi=static_cast<double>(j)-cogPhiLocal;
      row.moment_denominator+=energy;
      if (i!=3||j!=3)
      {
        row.moment_eta_numerator+=energy*deta*deta;
        row.moment_phi_numerator+=energy*dphi*dphi;
      }
      if (di<=1&&dj<=1)
      {
        row.moment33_denominator+=energy;
        if (i!=3||j!=3)
        {
          row.moment33_eta_numerator+=energy*deta*deta;
          row.moment33_phi_numerator+=energy*dphi*dphi;
        }
      }
    }
  }

  if (row.e33>0.0) row.e11_over_e33=row.e11/row.e33;
  if (row.e35>0.0) row.e32_over_e35=row.e32/row.e35;
  if (row.moment_denominator>0.0)
  {
    row.weta_cogx=row.moment_eta_numerator/row.moment_denominator;
    row.wphi_cogx=row.moment_phi_numerator/row.moment_denominator;
  }
  if (row.moment33_denominator>0.0)
  {
    row.weta33_cogx=row.moment33_eta_numerator/row.moment33_denominator;
    row.wphi33_cogx=row.moment33_phi_numerator/row.moment33_denominator;
  }
  row.finite_feature_state=
      std::isfinite(row.weta_cogx)&&std::isfinite(row.wphi_cogx)&&
      std::isfinite(row.weta33_cogx)&&std::isfinite(row.wphi33_cogx)&&
      std::isfinite(row.e11_over_e33)&&std::isfinite(row.e32_over_e35)&&
      std::all_of(nativeEt.begin(),nativeEt.end(),[](double item){return std::isfinite(item);});
  return row;
}

inline std::vector<float> modelFeatures(const ShowerFeatureViewRow& row,
                                        double photonEt,
                                        double vertexZ,
                                        double eta,
                                        double centrality,
                                        bool isAuAu)
{
  std::vector<float> result{
      static_cast<float>(photonEt),static_cast<float>(row.weta_cogx),
      static_cast<float>(row.wphi_cogx),static_cast<float>(vertexZ),
      static_cast<float>(eta),static_cast<float>(row.e11_over_e33),
      static_cast<float>(row.native_et1),static_cast<float>(row.native_et2),
      static_cast<float>(row.native_et3),static_cast<float>(row.native_et4),
      static_cast<float>(row.e32_over_e35)};
  if (isAuAu)
  {
    result.push_back(static_cast<float>(centrality));
    result.push_back(static_cast<float>(row.weta33_cogx));
    result.push_back(static_cast<float>(row.wphi33_cogx));
  }
  return result;
}
} // namespace RJShowerFactorialV1

#endif
