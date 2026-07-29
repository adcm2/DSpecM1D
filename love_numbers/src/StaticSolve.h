#ifndef DSPECM1D_LOVE_NUMBERS_STATIC_SOLVE_H
#define DSPECM1D_LOVE_NUMBERS_STATIC_SOLVE_H

#include <stdexcept>

#include <Eigen/SparseLU>

#include <DSpecM1D/LoveNumbers>

#include "LoveDofMap.h"
#include "RadialState.h"
#include "StaticForcing.h"
#include "StaticOperator.h"

namespace DSpecM1D::LoveNumbers::detail {

inline DegreeResult
solveStaticDegree(const EarthModels::ModelInput<double> &model,
                  RadialState &state, int degree) {
  if (degree < 0) {
    throw std::invalid_argument(
        "Static solve requires a non-negative degree.");
  }

  const LoveDofMap dofs(state.mesh(), degree);
  const Eigen::SparseMatrix<double> matrix =
      assembleStaticOperator(state, dofs, degree);
  const Eigen::MatrixXd right_hand_sides =
      assembleStaticRightHandSides(state, dofs, degree);

  Eigen::SparseLU<Eigen::SparseMatrix<double>> factorization;
  factorization.analyzePattern(matrix);
  factorization.factorize(matrix);
  if (factorization.info() != Eigen::Success) {
    throw std::runtime_error("Static operator factorization failed.");
  }

  const Eigen::MatrixXd solutions =
      factorization.solve(right_hand_sides);
  if (factorization.info() != Eigen::Success) {
    throw std::runtime_error("Static solve failed.");
  }

  const int surface_element = state.mesh().NE() - 1;
  const int surface_node = state.mesh().NN() - 1;
  const int surface_u =
      *dofs.index(Field::U, surface_element, surface_node);
  const auto surface_p =
      dofs.index(Field::P, surface_element, surface_node);
  const double density_norm = model.DensityNorm();
  const double length_norm = model.LengthNorm();
  const double time_norm = model.TimeNorm();
  const double potential_scale =
      length_norm / (density_norm * time_norm * time_norm);

  return DegreeResult{
      .degree = degree,
      .h_u = solutions(surface_u, 0) / density_norm,
      .k_u =
          surface_p ? solutions(*surface_p, 0) * potential_scale : 0.0,
      .h_phi = solutions(surface_u, 1) / density_norm,
      .k_phi =
          surface_p ? solutions(*surface_p, 1) * potential_scale : 0.0,
      .h_t =
          solutions(surface_u, 2) * time_norm * time_norm / length_norm,
      .k_t = surface_p ? solutions(*surface_p, 2) : 0.0,
  };
}

}   // namespace DSpecM1D::LoveNumbers::detail

#endif   // DSPECM1D_LOVE_NUMBERS_STATIC_SOLVE_H
