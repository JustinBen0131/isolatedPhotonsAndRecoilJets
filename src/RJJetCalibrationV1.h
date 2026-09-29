#ifndef RJ_JET_CALIBRATION_V1_H
#define RJ_JET_CALIBRATION_V1_H

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

// Method selection only, not a physics-approval or binary-compatibility gate.
// The qualification runner must separately bind matching headers/loaded DSO,
// payload bytes and collection-specific applicability before production.
namespace rj_jes_v1
{
inline constexpr const char* legacyCdbKey = "JES_Calib_Default";

template <class T, class = void>
struct HasExplicitMethodSelector : std::false_type {};

template <class T>
struct HasExplicitMethodSelector<T, std::void_t<decltype(
    std::declval<T&>().set_UseEMfracCalib(false))>> : std::true_type {};

template <class T>
void selectLegacyMethod(T& calibrator, const std::string& resolvedPayload,
                        const std::string& rawNode, const std::string& outputNode)
{
  if (resolvedPayload.empty())
    throw std::runtime_error("JES: resolved legacy payload is empty; unity fallback forbidden");
  if (rawNode.empty() || outputNode.empty() || rawNode == outputNode)
    throw std::runtime_error("JES: distinct pre-JES input and corrected output nodes required");
  if constexpr (HasExplicitMethodSelector<T>::value)
  {
    // In the newer JetCalib API set_CalibFile is EMfrac-only. Do not use it
    // to pretend that a legacy payload has been pinned.
    calibrator.set_UseEMfracCalib(false);
  }
  else
  {
    // Legacy-only headers cannot safely describe an EMfrac-default loaded
    // library. Refuse the mixed/unknown runtime instead of guessing its ABI.
    throw std::runtime_error("JES: explicit method-selector API unavailable; qualify matching JetCalib headers/library");
  }
}
}  // namespace rj_jes_v1

#endif
