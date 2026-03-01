// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TEST_SPEC_PRECONDITIONER_H
#define TEST_SPEC_PRECONDITIONER_H

#include <Eigen/IterativeLinearSolvers>
// #include "SolveBlock.h"

namespace Eigen {

template <typename Scalar_> class LUPD {
  typedef Scalar_ Scalar;
  typedef Matrix<Scalar, Dynamic, 1> Vector;

public:
  typedef typename Vector::StorageIndex StorageIndex;
  enum { ColsAtCompileTime = Dynamic, MaxColsAtCompileTime = Dynamic };

  LUPD() : m_isInitialized(false) {}

  template <typename MatType> explicit LUPD(const MatType &mat) {
    compute(mat);
  }

  constexpr Index rows() const noexcept { return matsize[0]; }
  constexpr Index cols() const noexcept { return matsize[1]; }

  template <typename MatType> LUPD &analyzePattern(const MatType &mat) {
    return *this;
  }

  template <typename MatType> LUPD &factorize(const MatType &mat) {
    // matsize[0] = mat.rows();
    // matsize[1] = mat.cols();
    // {
    //   Eigen::SparseMatrix<Scalar> smid;
    //   smid.resize(matsize[0], matsize[1]);
    //   smid.setIdentity();
    //   lu.compute(smid);
    // }
    m_isInitialized = true;
    return *this;
  }

  template <typename MatType> LUPD &compute(const MatType &mat) {
    // std::cout << "Hello 2\n";
    return analyzePattern(mat).factorize(mat);
  }

  // add matrix only invoked if we actually have a preconditioner. Otherwise it
  // will not have a preconditioner, with only the identity matrix used as the
  // preconditioner
  template <typename MatType> LUPD &addmatrix(const MatType &mat) {
    lu.compute(mat);
    return *this;
  }

  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs &b, Dest &x) const {
    x = lu.solve(b);
  }

  // solve call required by IterativeSolverBase
  template <typename Rhs>
  inline const Solve<LUPD, Rhs> solve(const MatrixBase<Rhs> &b) const {
    eigen_assert(m_isInitialized && "LUPD is not initialized.");
    eigen_assert(matsize[0] == b.rows() &&
                 "LUPD::solve(): invalid number of "
                 "rows of the right hand side matrix b");
    return Solve<LUPD, Rhs>(*this, b.derived());
  }

  ComputationInfo info() { return lu.info(); }

protected:
  // Vector m_invdiag;
  // Matrix<Scalar, Dynamic, Dynamic> m_val;
  std::vector<std::size_t> matsize{0, 0};
  bool m_isInitialized;
  SparseLU<SparseMatrix<Scalar>, COLAMDOrdering<int>> lu;
};

}   // namespace Eigen

#endif   // TEST_SPEC_PRECONDITIONER_H