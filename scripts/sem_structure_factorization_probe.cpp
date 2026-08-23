#include "config.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/SparseLU>
#include <lapacke.h>

#include <DSpecM1D/FrequencyTools>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/FullSpecMultiSem.h>
#include <DSpecM1D/src/LapackBandSolver.h>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>

// Benchmark-only exact static-condensation prototype.  It deliberately uses
// the assembled production matrix as its element-block source: element
// interiors are disjoint, while the assembled A_bb block is extracted once so
// shared endpoint contributions cannot be double counted.
namespace {
using C = std::complex<double>;
using Sparse = Eigen::SparseMatrix<C>;
using Dense = Eigen::Matrix<C, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point a, Clock::time_point b) {
  return std::chrono::duration<double>(b - a).count();
}

double infNorm(const Dense &x) {
  double result = 0;
  for (Eigen::Index row = 0; row < x.rows(); ++row) {
    double sum = 0;
    for (Eigen::Index col = 0; col < x.cols(); ++col)
      sum += std::abs(x(row, col));
    result = std::max(result, sum);
  }
  return result;
}

double infNorm(const Sparse &a) {
  Eigen::VectorXd rows = Eigen::VectorXd::Zero(a.rows());
  for (Eigen::Index col = 0; col < a.outerSize(); ++col)
    for (Sparse::InnerIterator it(a, col); it; ++it)
      rows(it.row()) += std::abs(it.value());
  return rows.size() ? rows.maxCoeff() : 0;
}

double residual(const Sparse &a, const Dense &x, const Dense &b) {
  return infNorm(a * x - b) /
         std::max(infNorm(a) * infNorm(x) + infNorm(b), 1.e-300);
}

double relative(const Dense &a, const Dense &b) {
  return infNorm(a - b) / std::max(infNorm(b), 1.e-300);
}

bool finite(const Dense &x) {
  return x.array().real().isFinite().all() &&
         x.array().imag().isFinite().all();
}

Dense rhs(Eigen::Index n) {
  Dense b(n, 3);
  for (Eigen::Index col = 0; col < b.cols(); ++col)
    for (Eigen::Index row = 0; row < b.rows(); ++row)
      b(row, col) = C(.35 + .0007 * row + .04 * col,
                      -.12 + .0003 * row - .02 * col);
  return b;
}

struct NaturalLU : Eigen::SparseLU<Sparse, Eigen::NaturalOrdering<int>> {
  NaturalLU() { this->setPivotThreshold(1.0); }
};

struct System {
  std::string family;
  std::string name;
  int degree = 0;
  int nq = 0;
  double maxstep = 0;
  double frequencyMhz = 0;
  Full1D::SEM sem;
  Sparse h, p, ha;
};

System makeSystem(const EarthModels::ModelInput<double, int> &model,
                  std::string family, std::string name, int degree, int nq,
                  double maxstep, double frequencyMhz) {
  System s;
  s.family = std::move(family);
  s.name = std::move(name);
  s.degree = degree;
  s.nq = nq;
  s.maxstep = maxstep;
  s.frequencyMhz = frequencyMhz;
  s.sem = Full1D::SEM(model, maxstep, nq, std::max(degree, 1));
  if (s.family == "radial") {
    s.h = s.sem.hR().cast<C>();
    s.p = s.sem.pR().cast<C>();
    s.ha = s.sem.hRa().cast<C>();
  } else if (s.family == "toroidal") {
    s.h = s.sem.hTk(degree).cast<C>();
    s.p = s.sem.pTk(degree).cast<C>();
    s.ha = s.sem.hTa(degree).cast<C>();
  } else {
    s.h = s.sem.hS(degree).cast<C>();
    s.p = s.sem.pS(degree).cast<C>();
    s.ha = s.sem.hSa(degree).cast<C>();
  }
  s.h.makeCompressed();
  s.p.makeCompressed();
  s.ha.makeCompressed();
  return s;
}

struct DynamicSystem {
  Sparse a;
  C wSquared;
  C chi;
};

DynamicSystem dynamicMatrix(const System &s,
                            const EarthModels::ModelInput<double, int> &model) {
  SpectraSolver::FreqFull grid(.1, 55.0, 1.0, 4.0, 1.0, .05, 0.0, 4.0, 1,
                               model.TimeNorm());
  SPARSESPEC::SpecConstants constants(grid.ep(), model.TREF());
  const double wreal = 2.0 * EIGEN_PI * s.frequencyMhz / 1000.0 *
                       model.TimeNorm();
  const C w = wreal + constants.ieps;
  const C chi = SPARSESPEC::attenFactor(wreal, constants.w0,
                                        constants.twodivpi, constants.myi);
  Sparse a = s.h - w * w * s.p;
  a += chi * s.ha;
  a.makeCompressed();
  return {std::move(a), w * w, chi};
}

struct Topology {
  std::vector<std::vector<Eigen::Index>> elementAll;
  std::vector<std::vector<Eigen::Index>> elementInterior;
  std::vector<std::vector<Eigen::Index>> elementBoundary;
  std::vector<Eigen::Index> boundary;
  std::vector<Eigen::Index> boundaryPosition;
  int firstElement = 0;
  int endElement = 0;
  int components = 0;
  int fluidSolidSplitDofs = 0;
};

Topology topology(const System &s) {
  Topology t;
  const int ne = s.sem.mesh().NE();
  t.firstElement = s.family == "toroidal" ? s.sem.el() : 0;
  t.endElement = s.family == "toroidal" ? s.sem.eu() : ne;
  t.components = s.family == "toroidal" ? 1 : (s.family == "radial" ? 2 : 3);
  std::set<Eigen::Index> retained;
  auto index = [&](int component, int element, int node) -> Eigen::Index {
    if (s.family == "radial") return s.sem.ltgR(component, element, node);
    if (s.family == "toroidal") return s.sem.ltgT(element, node);
    return s.sem.ltgS(component, element, node);
  };
  for (int element = t.firstElement; element < t.endElement; ++element) {
    std::vector<Eigen::Index> all, interior, boundary;
    for (int node = 0; node < s.nq; ++node)
      for (int component = 0; component < t.components; ++component) {
        const auto global = index(component, element, node);
        all.push_back(global);
        if (node == 0 || node == s.nq - 1) {
          boundary.push_back(global);
          retained.insert(global);
        } else {
          interior.push_back(global);
        }
      }
    t.elementAll.push_back(std::move(all));
    t.elementInterior.push_back(std::move(interior));
    t.elementBoundary.push_back(std::move(boundary));
  }
  t.boundary.assign(retained.begin(), retained.end());
  t.boundaryPosition.assign(static_cast<std::size_t>(s.h.rows()), -1);
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(t.boundary.size()); ++i)
    t.boundaryPosition[static_cast<std::size_t>(t.boundary[i])] = i;

  // A normal interface contributes one endpoint DOF per component.  Every
  // extra retained endpoint is the split tangential DOF at a fluid/solid
  // spheroidal interface.
  const int normal = t.components *
      (static_cast<int>(t.elementAll.size()) + 1);
  t.fluidSolidSplitDofs = static_cast<int>(t.boundary.size()) - normal;

  std::vector<int> owner(static_cast<std::size_t>(s.h.rows()), -1);
  for (std::size_t e = 0; e < t.elementInterior.size(); ++e)
    for (auto global : t.elementInterior[e]) {
      if (retained.count(global) || owner[static_cast<std::size_t>(global)] >= 0)
        throw std::runtime_error("interior DOF is not element-private");
      owner[static_cast<std::size_t>(global)] = static_cast<int>(e);
    }
  for (Eigen::Index col = 0; col < s.h.rows(); ++col) {
    const bool isBoundary = retained.count(col) != 0;
    if (!isBoundary && owner[static_cast<std::size_t>(col)] < 0)
      throw std::runtime_error("global DOF is absent from topology");
  }
  return t;
}

void verifyLocality(const Sparse &a, const Topology &t) {
  std::vector<int> owner(static_cast<std::size_t>(a.rows()), -1);
  for (std::size_t e = 0; e < t.elementInterior.size(); ++e)
    for (auto index : t.elementInterior[e])
      owner[static_cast<std::size_t>(index)] = static_cast<int>(e);
  for (Eigen::Index col = 0; col < a.outerSize(); ++col)
    for (Sparse::InnerIterator it(a, col); it; ++it) {
      const int rowOwner = owner[static_cast<std::size_t>(it.row())];
      const int colOwner = owner[static_cast<std::size_t>(it.col())];
      if (rowOwner >= 0 && colOwner >= 0 && rowOwner != colOwner &&
          std::abs(it.value()) != 0)
        throw std::runtime_error("nonlocal interior-interior matrix entry");
      if (rowOwner >= 0 && colOwner < 0) {
        bool found = false;
        for (auto boundary : t.elementBoundary[static_cast<std::size_t>(rowOwner)])
          found = found || boundary == it.col();
        if (!found && std::abs(it.value()) != 0)
          throw std::runtime_error("interior couples to nonlocal interface");
      }
      if (colOwner >= 0 && rowOwner < 0) {
        bool found = false;
        for (auto boundary : t.elementBoundary[static_cast<std::size_t>(colOwner)])
          found = found || boundary == it.row();
        if (!found && std::abs(it.value()) != 0)
          throw std::runtime_error("interface couples to nonlocal interior");
      }
      if (rowOwner < 0 && colOwner < 0 && it.row() != it.col() &&
          std::abs(it.value()) != 0) {
        bool found = false;
        for (const auto &elementBoundary : t.elementBoundary) {
          const bool hasRow = std::find(elementBoundary.begin(),
                                        elementBoundary.end(), it.row()) !=
                              elementBoundary.end();
          const bool hasCol = std::find(elementBoundary.begin(),
                                        elementBoundary.end(), it.col()) !=
                              elementBoundary.end();
          found = found || (hasRow && hasCol);
        }
        if (!found)
          throw std::runtime_error("nonlocal interface-interface matrix entry");
      }
    }
}

Dense extract(const Sparse &a, const std::vector<Eigen::Index> &rows,
              const std::vector<Eigen::Index> &cols) {
  Dense result(rows.size(), cols.size());
  for (Eigen::Index col = 0; col < result.cols(); ++col)
    for (Eigen::Index row = 0; row < result.rows(); ++row)
      result(row, col) = a.coeff(rows[static_cast<std::size_t>(row)],
                                 cols[static_cast<std::size_t>(col)]);
  return result;
}

Dense gather(const Dense &x, const std::vector<Eigen::Index> &rows) {
  Dense result(rows.size(), x.cols());
  for (Eigen::Index row = 0; row < result.rows(); ++row)
    result.row(row) = x.row(rows[static_cast<std::size_t>(row)]);
  return result;
}

void scatter(Dense &x, const std::vector<Eigen::Index> &rows, const Dense &v) {
  for (Eigen::Index row = 0; row < v.rows(); ++row)
    x.row(rows[static_cast<std::size_t>(row)]) = v.row(row);
}

struct LocalBlock {
  std::vector<Eigen::Index> interior;
  std::vector<Eigen::Index> boundaryGlobal;
  std::vector<Eigen::Index> boundaryReduced;
  Dense lu, aib, abi;
  std::vector<lapack_int> pivots;
};

struct LocalOperators {
  std::vector<Eigen::Index> interior;
  std::vector<Eigen::Index> boundaryGlobal;
  std::vector<Eigen::Index> boundaryReduced;
  Dense hii, pii, haii, hib, pib, haib, hbi, pbi, habi;
};

struct PreparedOperators {
  std::vector<LocalOperators> local;
  Dense hbb, pbb, habb;
  double extractionSeconds = 0;
};

PreparedOperators prepareOperators(const System &system, const Topology &top) {
  PreparedOperators prepared;
  const auto start = Clock::now();
  prepared.hbb = extract(system.h, top.boundary, top.boundary);
  prepared.pbb = extract(system.p, top.boundary, top.boundary);
  prepared.habb = extract(system.ha, top.boundary, top.boundary);
  for (std::size_t e = 0; e < top.elementInterior.size(); ++e) {
    LocalOperators local;
    local.interior = top.elementInterior[e];
    local.boundaryGlobal = top.elementBoundary[e];
    for (auto global : local.boundaryGlobal)
      local.boundaryReduced.push_back(
          top.boundaryPosition[static_cast<std::size_t>(global)]);
    local.hii = extract(system.h, local.interior, local.interior);
    local.pii = extract(system.p, local.interior, local.interior);
    local.haii = extract(system.ha, local.interior, local.interior);
    local.hib = extract(system.h, local.interior, local.boundaryGlobal);
    local.pib = extract(system.p, local.interior, local.boundaryGlobal);
    local.haib = extract(system.ha, local.interior, local.boundaryGlobal);
    local.hbi = extract(system.h, local.boundaryGlobal, local.interior);
    local.pbi = extract(system.p, local.boundaryGlobal, local.interior);
    local.habi = extract(system.ha, local.boundaryGlobal, local.interior);
    prepared.local.push_back(std::move(local));
  }
  prepared.extractionSeconds = seconds(start, Clock::now());
  return prepared;
}

lapack_int factor(Dense &a, std::vector<lapack_int> &pivots) {
  pivots.resize(static_cast<std::size_t>(std::min(a.rows(), a.cols())));
  return LAPACKE_zgetrf(LAPACK_COL_MAJOR, static_cast<lapack_int>(a.rows()),
                        static_cast<lapack_int>(a.cols()), a.data(),
                        static_cast<lapack_int>(a.rows()), pivots.data());
}

lapack_int solve(const Dense &lu, const std::vector<lapack_int> &pivots,
                 Dense &b) {
  return LAPACKE_zgetrs(LAPACK_COL_MAJOR, 'N',
                        static_cast<lapack_int>(lu.rows()),
                        static_cast<lapack_int>(b.cols()), lu.data(),
                        static_cast<lapack_int>(lu.rows()), pivots.data(),
                        b.data(), static_cast<lapack_int>(b.rows()));
}

struct Times {
  double dynamicPrep = 0, localFactor = 0, schur = 0, assembly = 0;
  double reducedPrep = 0;
  double interfaceFactor = 0, rhsCondense = 0, interfaceSolve = 0;
  double recovery = 0;
};

struct CondensedResult {
  Dense x;
  Times times;
  int info = 0;
  bool factorsFinite = true;
  Eigen::Index reducedKl = 0, reducedKu = 0;
  double localFlops = 0, schurFlops = 0, reducedFlops = 0;
};

CondensedResult condensedSolve(const PreparedOperators &operators,
                               C wSquared, C chi, Eigen::Index fullSize,
                               const Dense &b, const Topology &top) {
  CondensedResult result;
  std::vector<LocalBlock> blocks;
  blocks.reserve(top.elementInterior.size());

  auto start = Clock::now();
  Dense schur = operators.hbb - wSquared * operators.pbb +
                 chi * operators.habb;
  for (const auto &local : operators.local) {
    LocalBlock block;
    block.interior = local.interior;
    block.boundaryGlobal = local.boundaryGlobal;
    block.boundaryReduced = local.boundaryReduced;
    block.lu = local.hii - wSquared * local.pii + chi * local.haii;
    block.aib = local.hib - wSquared * local.pib + chi * local.haib;
    block.abi = local.hbi - wSquared * local.pbi + chi * local.habi;
    blocks.push_back(std::move(block));
  }
  auto end = Clock::now();
  result.times.dynamicPrep += seconds(start, end);

  for (auto &block : blocks) {
    const double ni = static_cast<double>(block.lu.rows());
    const double nb = static_cast<double>(block.aib.cols());
    start = Clock::now();
    result.info = factor(block.lu, block.pivots);
    end = Clock::now();
    result.times.localFactor += seconds(start, end);
    if (result.info != 0) return result;
    result.factorsFinite = result.factorsFinite && finite(block.lu);
    result.localFlops += 8.0 / 3.0 * ni * ni * ni;

    start = Clock::now();
    Dense solvedAib = block.aib;
    result.info = solve(block.lu, block.pivots, solvedAib);
    if (result.info != 0) return result;
    Dense update = block.abi * solvedAib;
    end = Clock::now();
    result.times.schur += seconds(start, end);
    result.schurFlops += 8.0 * ni * ni * nb + 8.0 * nb * ni * nb;

    start = Clock::now();
    for (Eigen::Index col = 0; col < update.cols(); ++col)
      for (Eigen::Index row = 0; row < update.rows(); ++row)
        schur(block.boundaryReduced[static_cast<std::size_t>(row)],
              block.boundaryReduced[static_cast<std::size_t>(col)]) -=
            update(row, col);
    end = Clock::now();
    result.times.assembly += seconds(start, end);
  }

  for (Eigen::Index col = 0; col < schur.cols(); ++col)
    for (Eigen::Index row = 0; row < schur.rows(); ++row)
      if (schur(row, col) != C{}) {
        result.reducedKl = std::max(result.reducedKl, row - col);
        result.reducedKu = std::max(result.reducedKu, col - row);
      }
  SPARSESPEC::detail::LapackBandSolver interfaceSolver;
  start = Clock::now();
  auto &band = interfaceSolver.bandWorkspace(schur.rows(), result.reducedKl,
                                              result.reducedKu);
  band.data.setZero();
  for (Eigen::Index col = 0; col < schur.cols(); ++col)
    for (Eigen::Index row = std::max<Eigen::Index>(0, col - result.reducedKu);
         row <= std::min<Eigen::Index>(schur.rows() - 1,
                                      col + result.reducedKl); ++row)
      band.at(row, col) = schur(row, col);
  end = Clock::now();
  result.times.reducedPrep += seconds(start, end);
  start = Clock::now();
  interfaceSolver.factorize();
  end = Clock::now();
  result.times.interfaceFactor += seconds(start, end);
  result.info = interfaceSolver.info();
  if (result.info != 0) return result;
  result.factorsFinite = result.factorsFinite &&
      finite(interfaceSolver.bandWorkspace().data);
  const double nr = static_cast<double>(schur.rows());
  result.reducedFlops = 8.0 * nr *
      static_cast<double>(result.reducedKl + result.reducedKu + 1) *
      static_cast<double>(result.reducedKl + 1);

  start = Clock::now();
  Dense condensed = gather(b, top.boundary);
  for (const auto &block : blocks) {
    Dense yi = gather(b, block.interior);
    result.info = solve(block.lu, block.pivots, yi);
    if (result.info != 0) return result;
    const Dense contribution = block.abi * yi;
    for (Eigen::Index row = 0; row < contribution.rows(); ++row)
      condensed.row(block.boundaryReduced[static_cast<std::size_t>(row)]) -=
          contribution.row(row);
  }
  end = Clock::now();
  result.times.rhsCondense += seconds(start, end);

  start = Clock::now();
  Dense xb = interfaceSolver.solve(condensed);
  end = Clock::now();
  result.times.interfaceSolve += seconds(start, end);
  result.info = interfaceSolver.info();
  if (result.info != 0) return result;

  start = Clock::now();
  result.x = Dense::Zero(fullSize, b.cols());
  scatter(result.x, top.boundary, xb);
  for (const auto &block : blocks) {
    Dense localBoundary(block.boundaryGlobal.size(), b.cols());
    for (Eigen::Index row = 0; row < localBoundary.rows(); ++row)
      localBoundary.row(row) =
          xb.row(block.boundaryReduced[static_cast<std::size_t>(row)]);
    Dense interiorRhs = gather(b, block.interior) - block.aib * localBoundary;
    result.info = solve(block.lu, block.pivots, interiorRhs);
    if (result.info != 0) return result;
    scatter(result.x, block.interior, interiorRhs);
  }
  end = Clock::now();
  result.times.recovery += seconds(start, end);
  return result;
}

struct ReferenceResult {
  Dense x;
  double factorSeconds = 0, solveSeconds = 0;
  int info = 0;
};

ReferenceResult eigenReference(const Sparse &a, const Dense &b) {
  ReferenceResult r;
  NaturalLU solver;
  solver.analyzePattern(a);
  auto start = Clock::now();
  solver.factorize(a);
  auto end = Clock::now();
  r.factorSeconds = seconds(start, end);
  r.info = solver.info();
  if (r.info != Eigen::Success) return r;
  start = Clock::now();
  r.x = solver.solve(b);
  end = Clock::now();
  r.solveSeconds = seconds(start, end);
  r.info = solver.info();
  return r;
}

ReferenceResult lapackReference(const Sparse &a, const Dense &b) {
  ReferenceResult r;
  const auto packed = SPARSESPEC::detail::packLapackBand(a);
  SPARSESPEC::detail::LapackBandSolver solver;
  auto &workspace = solver.bandWorkspace(packed.n, packed.kl, packed.ku);
  workspace.data = packed.data;
  auto start = Clock::now();
  solver.factorize();
  auto end = Clock::now();
  r.factorSeconds = seconds(start, end);
  r.info = solver.info();
  if (r.info != 0) return r;
  start = Clock::now();
  r.x = solver.solve(b);
  end = Clock::now();
  r.solveSeconds = seconds(start, end);
  r.info = solver.info();
  return r;
}

void add(Times &sum, const Times &value) {
  sum.dynamicPrep += value.dynamicPrep;
  sum.localFactor += value.localFactor;
  sum.schur += value.schur;
  sum.assembly += value.assembly;
  sum.reducedPrep += value.reducedPrep;
  sum.interfaceFactor += value.interfaceFactor;
  sum.rhsCondense += value.rhsCondense;
  sum.interfaceSolve += value.interfaceSolve;
  sum.recovery += value.recovery;
}

void divide(Times &value, double n) {
  value.dynamicPrep /= n; value.localFactor /= n; value.schur /= n;
  value.assembly /= n; value.reducedPrep /= n; value.interfaceFactor /= n;
  value.rhsCondense /= n; value.interfaceSolve /= n; value.recovery /= n;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: sem_structure_factorization_probe MODEL\n";
    return 2;
  }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    std::vector<System> systems;
    systems.push_back(makeSystem(model, "radial", "radial_0S0", 0, 6,
                                 .05, .8146635));
    systems.push_back(makeSystem(model, "toroidal", "toroidal_0T12", 12, 6,
                                 .05, 1.881496));
    systems.push_back(makeSystem(model, "spheroidal", "spheroidal_0S12", 12,
                                 6, .05, 2.003070));
    systems.push_back(makeSystem(model, "fluid-solid", "fluid_solid_0S40", 40,
                                 6, .02, 4.761725));
    systems.push_back(makeSystem(model, "large-spheroidal", "large_0S40", 40,
                                 6, .01, 4.761725));

    std::cout << std::setprecision(17)
      << "family\tname\tdegree\tfrequency_mhz\tnq\tmaxstep\telements\tn\tnnz"
      << "\tkl\tku\tcomponents\ttotal_interior\tinterface_dofs\tfluid_solid_split_dofs"
      << "\treduced_n\treduced_kl\treduced_ku\tlocal_interior_min\tlocal_interior_max"
      << "\tlocal_lu_flops\tschur_flops\treduced_factor_flops\trepeats"
      << "\tone_time_base_extraction_s\tlocal_dynamic_prep_s"
      << "\tlocal_factor_s\tschur_s\treduced_assembly_s"
      << "\treduced_band_prep_s\tinterface_factor_s"
      << "\trhs_condense_s\tinterface_solve_s\trecovery_s\tcondensed_setup_intrinsic_s"
      << "\tcondensed_setup_with_scaffolding_s\tcondensed_solve_recovery_s"
      << "\teigen_factor_s\teigen_solve_s\tlapack_factor_s\tlapack_solve_s"
      << "\tsetup_speedup_vs_lapack\ttotal_speedup_vs_lapack"
      << "\tcondensed_info\teigen_info\tlapack_info\tfactors_finite\tfinite"
      << "\tcondensed_residual\teigen_residual\tlapack_residual"
      << "\tcondensed_vs_eigen\tcondensed_vs_lapack\teigen_vs_lapack\n";

    constexpr int repeats = 3;
    volatile double sink = 0;
    for (const auto &system : systems) {
      const auto dynamic = dynamicMatrix(system, model);
      const Sparse &a = dynamic.a;
      const Dense b = rhs(a.rows());
      const auto top = topology(system);
      verifyLocality(a, top);
      const auto operators = prepareOperators(system, top);
      const auto bandwidth = SPARSESPEC::detail::lapackBandBandwidth(a);

      // Exactly one warm-up of each solver.
      auto warmCondensed = condensedSolve(operators, dynamic.wSquared,
                                           dynamic.chi, a.rows(), b, top);
      auto warmEigen = eigenReference(a, b);
      auto warmLapack = lapackReference(a, b);
      if (warmCondensed.x.size()) sink += warmCondensed.x(0, 0).real();
      if (warmEigen.x.size()) sink += warmEigen.x(0, 0).real();
      if (warmLapack.x.size()) sink += warmLapack.x(0, 0).real();

      Times condensedTimes;
      CondensedResult condensed;
      ReferenceResult eigen, lapack;
      double eigenFactor = 0, eigenSolve = 0, lapackFactor = 0, lapackSolve = 0;
      for (int repeat = 0; repeat < repeats; ++repeat) {
        condensed = condensedSolve(operators, dynamic.wSquared, dynamic.chi,
                                    a.rows(), b, top);
        add(condensedTimes, condensed.times);
        eigen = eigenReference(a, b);
        eigenFactor += eigen.factorSeconds;
        eigenSolve += eigen.solveSeconds;
        lapack = lapackReference(a, b);
        lapackFactor += lapack.factorSeconds;
        lapackSolve += lapack.solveSeconds;
      }
      divide(condensedTimes, repeats);
      eigenFactor /= repeats; eigenSolve /= repeats;
      lapackFactor /= repeats; lapackSolve /= repeats;

      const double intrinsicSetup = condensedTimes.dynamicPrep +
          condensedTimes.localFactor +
          condensedTimes.schur + condensedTimes.assembly +
          condensedTimes.interfaceFactor;
      const double withScaffolding = intrinsicSetup + condensedTimes.reducedPrep;
      const double condensedSolveTime = condensedTimes.rhsCondense +
          condensedTimes.interfaceSolve + condensedTimes.recovery;
      const double lapackTotal = lapackFactor + lapackSolve;
      const bool allFinite = finite(condensed.x) && finite(eigen.x) &&
                             finite(lapack.x);
      std::size_t totalInterior = 0;
      Eigen::Index localMin = std::numeric_limits<Eigen::Index>::max();
      Eigen::Index localMax = 0;
      for (const auto &local : top.elementInterior) {
        totalInterior += local.size();
        localMin = std::min(localMin, static_cast<Eigen::Index>(local.size()));
        localMax = std::max(localMax, static_cast<Eigen::Index>(local.size()));
      }
      std::cout
        << system.family << '\t' << system.name << '\t' << system.degree << '\t'
        << system.frequencyMhz << '\t' << system.nq << '\t' << system.maxstep << '\t'
        << top.elementAll.size() << '\t' << a.rows() << '\t' << a.nonZeros() << '\t'
        << bandwidth.first << '\t' << bandwidth.second << '\t' << top.components << '\t'
        << totalInterior << '\t' << top.boundary.size() << '\t'
        << top.fluidSolidSplitDofs << '\t' << top.boundary.size() << '\t'
        << condensed.reducedKl << '\t' << condensed.reducedKu << '\t'
        << localMin << '\t' << localMax << '\t' << condensed.localFlops << '\t'
        << condensed.schurFlops << '\t' << condensed.reducedFlops << '\t'
        << repeats << '\t' << operators.extractionSeconds << '\t'
        << condensedTimes.dynamicPrep << '\t' << condensedTimes.localFactor << '\t'
        << condensedTimes.schur << '\t'
        << condensedTimes.assembly << '\t' << condensedTimes.reducedPrep << '\t'
        << condensedTimes.interfaceFactor << '\t'
        << condensedTimes.rhsCondense << '\t' << condensedTimes.interfaceSolve << '\t'
        << condensedTimes.recovery << '\t' << intrinsicSetup << '\t'
        << withScaffolding << '\t' << condensedSolveTime << '\t'
        << eigenFactor << '\t' << eigenSolve << '\t' << lapackFactor << '\t'
        << lapackSolve << '\t' << lapackFactor / intrinsicSetup << '\t'
        << lapackTotal / (intrinsicSetup + condensedSolveTime) << '\t'
        << condensed.info << '\t' << eigen.info << '\t' << lapack.info << '\t'
        << condensed.factorsFinite << '\t' << allFinite << '\t'
        << residual(a, condensed.x, b) << '\t'
        << residual(a, eigen.x, b) << '\t' << residual(a, lapack.x, b) << '\t'
        << relative(condensed.x, eigen.x) << '\t'
        << relative(condensed.x, lapack.x) << '\t'
        << relative(eigen.x, lapack.x) << '\n';
    }
    std::cerr << "sink\t" << sink << '\n';
  } catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
