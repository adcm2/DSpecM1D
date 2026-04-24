#ifndef DSPECM1D_MODEL_INFO_CONCEPTS_H
#define DSPECM1D_MODEL_INFO_CONCEPTS_H

#include <concepts>

namespace ModelInfoConcepts {

template <typename T>
concept HasBasicNormalisationInformation = requires(T t) {
  { t.LengthNorm() } -> std::convertible_to<double>;
  { t.MassNorm() } -> std::convertible_to<double>;
  { t.TimeNorm() } -> std::convertible_to<double>;
};

template <typename Model, typename INT, typename FLOAT>
concept BasicSphericalGeometryModel = requires(Model model, INT i, FLOAT r) {
  requires std::integral<INT>;
  requires std::floating_point<FLOAT>;
  requires HasBasicNormalisationInformation<Model>;

  { model.NumberOfLayers() } -> std::convertible_to<INT>;
  { model.LowerRadius(i) } -> std::convertible_to<FLOAT>;
  { model.UpperRadius(i) } -> std::convertible_to<FLOAT>;
  { model.OuterRadius() } -> std::convertible_to<FLOAT>;
};

template <typename Model, typename INT, typename FLOAT>
concept BasicSphericalDensityModel = requires(Model model, INT i, FLOAT r) {
  requires BasicSphericalGeometryModel<Model, INT, FLOAT>;

  { model.Density(i) } -> std::regular_invocable<FLOAT>;
  { model.Density(i)(r) } -> std::convertible_to<FLOAT>;
};

}  // namespace ModelInfoConcepts

#endif
