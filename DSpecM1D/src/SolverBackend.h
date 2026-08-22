#ifndef DSPECM1D_SOLVER_BACKEND_H
#define DSPECM1D_SOLVER_BACKEND_H

namespace SPARSESPEC {

/// Selects the numerical backend used by the preferred adaptive solver.
enum class SolverBackend {
  EigenSparseLU,
  LapackBandLU
};

} // namespace SPARSESPEC

#endif // DSPECM1D_SOLVER_BACKEND_H
