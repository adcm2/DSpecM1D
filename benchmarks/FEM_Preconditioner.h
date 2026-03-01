// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef FEM_PRECONDITIONER_H
#define FEM_PRECONDITIONER_H

#include <Eigen/IterativeLinearSolvers>
// #include "SolveBlock.h"

namespace Eigen {

template <typename Scalar_> class FEMPreconditioner {
   typedef Scalar_ Scalar;
   typedef Matrix<Scalar, Dynamic, 1> Vector;

 public:
   typedef typename Vector::StorageIndex StorageIndex;
   enum { ColsAtCompileTime = Dynamic, MaxColsAtCompileTime = Dynamic };

   FEMPreconditioner() : m_isInitialized(false) {}

   template <typename MatType>
   explicit FEMPreconditioner(const MatType &mat) {
      compute(mat);
   }

   constexpr Index rows() const noexcept { return matsize[0]; }
   constexpr Index cols() const noexcept { return matsize[1]; }

   template <typename MatType>
   FEMPreconditioner &analyzePattern(const MatType &mat) {
      return *this;
   }

   template <typename MatType>
   FEMPreconditioner &factorize(const MatType &mat) {
      matsize[0] = mat.rows();
      matsize[1] = mat.cols();
      if (!m_isInitialized){
         Eigen::SparseMatrix<Scalar> smid;
         smid.resize(matsize[0], matsize[1]);
         smid.setIdentity();
         lu_solver.compute(smid);
      }
      m_isInitialized = true;
      return *this;
   }

   template <typename MatType>
   FEMPreconditioner &compute(const MatType &mat) {
      // std::cout << "Hello 2\n";
      return analyzePattern(mat).factorize(mat);
   }

   // add matrix only invoked if we actually have a preconditioner. Otherwise it
   // will not have a preconditioner, with only the identity matrix used as the
   // preconditioner
   template <typename MatType>
   FEMPreconditioner &addmatrix(const MatType &mat) {
      lu_solver.compute(mat);
      matsize[0] = mat.rows();
      matsize[1] = mat.cols();
      m_isInitialized = true;
      return *this;
   }

   template <typename Rhs, typename Dest>
   void _solve_impl(const Rhs &b, Dest &x) const {
      x = lu_solver.solve(b);
   }

   // solve call required by IterativeSolverBase
   template <typename Rhs>
   inline const Solve<FEMPreconditioner, Rhs>
   solve(const MatrixBase<Rhs> &b) const {
      eigen_assert(m_isInitialized &&
                   "FEMPreconditioner is not initialized.");
      eigen_assert(
          matsize[0] == b.rows() &&
          "FEMPreconditioner::solve(): invalid number of "
          "rows of the right hand side matrix b");
      return Solve<FEMPreconditioner, Rhs>(*this, b.derived());
   }

   ComputationInfo info() { return lu_solver.info(); }

 protected:
   // Vector m_invdiag;
   // Matrix<Scalar, Dynamic, Dynamic> m_val;
   std::vector<std::size_t> matsize{0, 0};
   bool m_isInitialized;
   SparseLU<SparseMatrix<Scalar>, Eigen::COLAMDOrdering<int>> lu_solver;
};

}   // namespace Eigen

#endif   // SPHERICAL_GEOMETRY_PRECONDITIONER_H