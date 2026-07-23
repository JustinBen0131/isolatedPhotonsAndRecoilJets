#include "../../../src/RJShowerFactorialV1.h"

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using namespace RJReplayFoundationV1;
using namespace RJShowerFactorialV1;

void require(bool condition,const std::string& message)
{
  if(!condition)throw std::runtime_error(message);
}

ShowerCellRow cell(const Identity128& candidate,int eta,int phi,double calibrated,
                   bool owned,double rawValue,int gridMask)
{
  ShowerCellRow row;
  row.candidate_id=candidate;
  row.tower_eta_index=eta;
  row.tower_phi_index=(phi+256)%256;
  row.local_eta_index=eta-49;
  int localPhi=row.tower_phi_index;
  while(localPhi>127)localPhi-=256;
  row.local_phi_index=localPhi;
  row.tower_key=static_cast<std::uint64_t>(eta*256+row.tower_phi_index+1);
  row.calibrated_energy=calibrated;
  row.rawcluster_owned=owned?1:0;
  row.rawcluster_value_present=owned?1:0;
  row.rawcluster_map_value=owned?rawValue:std::numeric_limits<double>::quiet_NaN();
  row.is_good=1;
  row.is_zero=calibrated==0.0;
  row.is_negative=calibrated<0.0;
  row.is_nonfinite=!std::isfinite(calibrated);
  row.denominator_membership=std::isfinite(calibrated)&&calibrated>0.0;
  row.floor0_membership=row.denominator_membership;
  row.floor70_membership=std::isfinite(calibrated)&&calibrated>0.070;
  row.grid_membership_bitmask=gridMask;
  return row;
}
}

int main()
{
  using namespace RJReplayFoundationV1;
  using namespace RJShowerFactorialV1;
  try
  {
    const Identity128 candidate=makeIdentity("THE134|synthetic|candidate");
    // H70 is centered across the CEMC phi seam from H0.  The retained cells
    // are an absolute-key union, so both definitions remain exactly replayable.
    const double h0Eta=48.25,h0Phi=255.75;
    const double h70Eta=49.25,h70Phi=0.25;
    std::vector<ShowerCellRow> cells;
    std::set<std::pair<int,int>> seen;
    auto appendGrid=[&](int centerEta,int centerPhi,int mask)
    {
      for(int de=-3;de<=3;++de)for(int dp=-3;dp<=3;++dp)
      {
        const int eta=centerEta+de,phi=(centerPhi+dp+256)%256;
        if(eta<0||eta>=96)continue;
        if(!seen.insert({eta,phi}).second)
        {
          for(auto& row:cells)if(row.tower_eta_index==eta&&row.tower_phi_index==phi)
            row.grid_membership_bitmask|=mask;
          continue;
        }
        cells.push_back(cell(candidate,eta,phi,0.0,false,0.0,mask));
      }
    };
    appendGrid(48,255,1);
    appendGrid(49,0,2);

    auto setCell=[&](int eta,int phi,double calibrated,bool owned,double raw)
    {
      phi=(phi+256)%256;
      for(auto& row:cells)if(row.tower_eta_index==eta&&row.tower_phi_index==phi)
      {
        const int mask=row.grid_membership_bitmask;
        row=cell(candidate,eta,phi,calibrated,owned,raw,mask);
        return;
      }
      throw std::runtime_error("synthetic tower is outside retained union");
    };
    setCell(49,0,2.0,true,1.7);       // H70 center, owned.
    setCell(50,0,0.5,true,0.3);       // owned, changes raw-map diagnostic.
    setCell(49,1,1.0,false,0.0);      // full-grid only, separates H/G/O.
    setCell(48,255,0.05,true,0.04);   // retained only by zero-floor views.

    const std::array<double,4> native0{{1.1,1.2,1.3,1.4}};
    const std::array<double,4> native70{{2.1,2.2,2.3,2.4}};
    std::map<std::string,ShowerFeatureViewRow> views;
    std::set<std::string> hashes;
    for(const auto& def:definitions())
    {
      const bool zeroFloor=def.floor_gev==0.0;
      auto view=buildView(candidate,def,cells,
                          zeroFloor?h0Eta:h70Eta,
                          zeroFloor?h0Phi:h70Phi,
                          zeroFloor?native0:native70);
      view.ordered_features=modelFeatures(view,20.0,3.0,0.1,30.0,true);
      require(view.ordered_features.size()==14,"AuAu feature order is not 14");
      require(hashes.insert(view.semantic_sha256).second,"definition semantic hashes are not unique");
      views.emplace(def.name,std::move(view));
    }

    require(views.size()==7,"factorial definition count is not seven");
    require(views.at("H70").center_phi_index==0,"H70 phi seam center was not normalized");
    require(views.at("H0").center_phi_index==255,"H0 center was not retained independently");
    require(views.at("H70").e33>views.at("O70").e33,"full-grid sums do not differ from ownership-only sums");
    require(views.at("H70").moment_denominator<views.at("G70").moment_denominator,
            "ownership-masked moments do not differ from full-grid moments");
    require(views.at("H0").moment_denominator>views.at("H70").moment_denominator,
            "zero-floor owned population is not represented");
    require(std::fabs(views.at("R70").moment_denominator-views.at("O70").moment_denominator)>1.0e-12,
            "raw-map and calibrated ownership views are indistinguishable");
    require(modelFeatures(views.at("H70"),20.0,3.0,0.1,30.0,false).size()==11,
            "pp feature order is not 11");
    const auto canonicalAuAuFeatures=modelFeatures(
        views.at("H70"),20.0,3.0,0.1,30.0,true);
    const std::array<float,14> expectedAuAuFeatures{{
        20.0F,
        static_cast<float>(views.at("H70").weta_cogx),
        static_cast<float>(views.at("H70").wphi_cogx),
        static_cast<float>(views.at("H70").weta33_cogx),
        static_cast<float>(views.at("H70").wphi33_cogx),
        3.0F,0.1F,
        static_cast<float>(views.at("H70").e11_over_e33),
        static_cast<float>(views.at("H70").native_et1),
        static_cast<float>(views.at("H70").native_et2),
        static_cast<float>(views.at("H70").native_et3),
        static_cast<float>(views.at("H70").native_et4),
        static_cast<float>(views.at("H70").e32_over_e35),
        30.0F}};
    require(canonicalAuAuFeatures.size()==expectedAuAuFeatures.size(),
            "AuAu feature order is not 14");
    require(std::equal(canonicalAuAuFeatures.begin(),canonicalAuAuFeatures.end(),
                       expectedAuAuFeatures.begin()),
            "AuAu canonical feature order drifted: 3x3 widths must follow cog widths and centrality must be last");
    for(const auto& item:views)
      require(item.second.finite_feature_state==1,"nonfinite factorial view: "+item.first);

    // Regression from the first real Photon20 canary mismatch.  The
    // authoritative PhotonClusterBuilder accumulates these terms in float32
    // row-major order.  Double accumulation changes wphi_cogx by 7.99e-6 and
    // is not an exact model-input replay even though it is numerically close.
    const Identity128 precisionCandidate=makeIdentity("THE134|float32-regression|candidate");
    std::vector<ShowerCellRow> precisionCells;
    for(int de=-3;de<=3;++de)for(int dp=-3;dp<=3;++dp)
      precisionCells.push_back(cell(precisionCandidate,61+de,64+dp,0.0,false,0.0,1));
    auto setPrecisionCell=[&](int eta,int phi,float energy)
    {
      for(auto& row:precisionCells)if(row.tower_eta_index==eta&&row.tower_phi_index==phi)
      {
        row=cell(precisionCandidate,eta,phi,energy,true,energy,1);
        return;
      }
      throw std::runtime_error("float32 regression tower is outside retained grid");
    };
    setPrecisionCell(60,63,0.09513379633426666F);
    setPrecisionCell(61,62,0.6544426083564758F);
    setPrecisionCell(61,63,0.8515732884407043F);
    setPrecisionCell(61,64,3.8511412143707275F);
    setPrecisionCell(61,65,0.07497962564229965F);
    setPrecisionCell(62,62,8.815872192382812F);
    setPrecisionCell(62,63,3.3401453495025635F);
    setPrecisionCell(62,64,2.131049633026123F);
    setPrecisionCell(63,63,0.10414200276136398F);
    const auto precisionView=buildView(
        precisionCandidate,definition("H70"),precisionCells,
        static_cast<double>(61.83404541015625F),
        static_cast<double>(64.3939323425293F),native70);
    require(static_cast<float>(precisionView.wphi_cogx)==1.8874231576919556F,
            "H70 wphi_cogx does not preserve PhotonClusterBuilder float32 arithmetic");

    std::cout<<"THE134_SHOWER_FACTORIAL_CANARY_PASS"
             <<" cells="<<cells.size()<<" views="<<views.size()
             <<" h70_hash="<<views.at("H70").semantic_sha256<<std::endl;
    return 0;
  }
  catch(const std::exception& error)
  {
    std::cerr<<"THE134_SHOWER_FACTORIAL_CANARY_FAIL: "<<error.what()<<std::endl;
    return 1;
  }
}
