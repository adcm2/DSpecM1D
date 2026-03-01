#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Spectra/SymEigsShiftSolver.h>
// <Spectra/MatOp/DenseSymShiftSolve.h> is implicitly included
#include <iostream>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <random>

using namespace Spectra;

int
main() {
  // A size-10 diagonal matrix with elements 1, 2, ..., 10
  int mlen = 10;
  int N = 8;
  if (N > mlen - 1) {
    N = mlen - 1;
  }
  long int maxn = 2 * N;
  if (maxn > mlen) {
    maxn = mlen;
  }

  // matrix m
  Eigen::MatrixXd MM = Eigen::MatrixXd::Zero(mlen, mlen);
  for (int i = 0; i < MM.rows(); i++)
    MM(i, i) = i + 1;

  // Construct matrix operation object using the wrapper class
  DenseSymShiftSolve<double> op_dense(MM);

  // Construct eigen solver object with shift 0
  // This will find eigenvalues that are closest to 0#
  double sigshift = 0.0;
  std::cout << "Enter shift:\n";
  std::cin >> sigshift;
  SymEigsShiftSolver<DenseSymShiftSolve<double>> eigs(op_dense, N, maxn,
                                                      sigshift);

  eigs.init();
  /*
  eigs.compute(SortRule::LargestMagn);
  if (eigs.info() == CompInfo::Successful) {
    Eigen::VectorXd evalues = eigs.eigenvalues();
    // Will get (3.0, 2.0, 1.0)
    std::cout << "Eigenvalues found:\n" << evalues << std::endl;
  }
  */

  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////

  // declaring necessary objects
  Eigen::SparseMatrix<double> m_ke(mlen, mlen), m_in(mlen, mlen),
      m_se(mlen, mlen);
  using T = Eigen::Triplet<double>;
  std::vector<T> tpl_ke, tpl_in, tpl_seig;
  std::vector<double> vec_nn, vec_lm;
  double lower_bound = 0.001;
  double upper_bound = 1.0;
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  std::default_random_engine re;
  for (int idx = 0; idx < mlen; ++idx) {
    vec_nn.push_back(unif(re));
    vec_lm.push_back(1.0 / vec_nn[idx]);
  }
  for (int idx = 0; idx < mlen; ++idx) {
    // vec_nn.push_back(unif(re));
    vec_lm.push_back(1.0 / vec_nn[idx]);
  }
  // for (auto i : vec_nn) {
  //   std::cout << i << "\n";
  // }
  // std::cout << std::endl;

  for (int i = 0; i < mlen; i++) {
    for (int j = 0; j < i + 1; j++) {
      double curval = unif(re);
      tpl_ke.push_back(T(i, j, curval));
      tpl_ke.push_back(T(j, i, curval));
      tpl_seig.push_back(T(i, j, curval * vec_lm[i]));
      tpl_seig.push_back(T(j, i, curval * vec_lm[j]));
    }

    tpl_in.push_back(T(i, i, vec_nn[i]));
  }

  m_ke.setFromTriplets(tpl_ke.begin(), tpl_ke.end());
  m_in.setFromTriplets(tpl_in.begin(), tpl_in.end());
  m_se.setFromTriplets(tpl_seig.begin(), tpl_seig.end());
  m_ke.makeCompressed();
  m_in.makeCompressed();
  m_se.makeCompressed();

  SparseGenMatProd<double> op_seig(m_se);

  GenEigsSolver<SparseGenMatProd<double>> eigs2(op_seig, N, maxn);
  eigs2.init();
  int n2 = eigs2.compute(SortRule::SmallestMagn);
  if (eigs2.info() == CompInfo::Successful) {
    // std::cout << "Successful\n";
    std::cout << "\nFrom SparseGenMatProd: \n";
    std::cout << eigs2.eigenvalues() << "\n\n";

  } else {
    std::cout << "Unsuccessful non-generalised\n";
  }
  //   std::cout << m_ke << "\n";
  //   std::cout << m_in << "\n";
  using OpType = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  using BOpType = SparseSymMatProd<double>;
  OpType op(m_ke, m_in);
  BOpType op_ke(m_ke), Bop(m_in);

  /////////////////////////////////////////////////////////////////////
  // solving eigenproblem
  //  initiate eigensolver
  SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> eig_gen(
      op, Bop, N, maxn, sigshift);
  eig_gen.init();
  int nconv_gen = eig_gen.compute(SortRule::LargestMagn);
  Eigen::VectorXcd evalues_gen;
  if (eig_gen.info() == CompInfo::Successful) {
    // std::cout << "Eigenvalues from generalised: \n"
    //           << eig_gen.eigenvalues() * _freq_norm * _freq_norm << "\n\n";
    evalues_gen = eig_gen.eigenvalues();
    std::cout << "\nFrom SymGEigsShiftSolver: \n";
    std::cout << evalues_gen << "\n\n";
  } else {
    std::cout << "Unsuccessful generalised\n";
  }

  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  /*
  const int n = 100;
  // Define the A matrix
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(n, n);
  Eigen::MatrixXd A = M + M.transpose();

  // Define the B matrix, a tridiagonal matrix with 2 on the diagonal
  // and 1 on the subdiagonals
  Eigen::SparseMatrix<double> B(n, n);
  B.reserve(Eigen::VectorXi::Constant(n, 3));
  for (int i = 0; i < n; i++) {
    B.insert(i, i) = 2.0;
    if (i > 0)
      B.insert(i - 1, i) = 1.0;
    if (i < n - 1)
      B.insert(i + 1, i) = 1.0;
  }
  Eigen::SparseMatrix<double> C = A.sparseView();

  // Construct matrix operation objects using the wrapper classes
  // A is dense, B is sparse
  using OpType = SymShiftInvert<double, Eigen::Dense, Eigen::Sparse>;
  using BOpType = SparseSymMatProd<double>;
  using OpType2 = SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
  using BOpType2 = SparseSymMatProd<double>;
  OpType op(A, B);
  BOpType Bop(B);
  OpType2 op2(A, B);
  BOpType2 Bop2(B);

  // Construct generalized eigen solver object, seeking three generalized
  // eigenvalues that are closest to zero. This is equivalent to specifying
  // a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule
  SymGEigsShiftSolver<OpType, BOpType, GEigsMode::ShiftInvert> geigs(op, Bop, 3,
                                                                     6, 0.0);
  SymGEigsShiftSolver<OpType2, BOpType2, GEigsMode::ShiftInvert> geigs2(
      op2, Bop2, 3, 6, 0.0);

  // Initialize and compute
  geigs.init();
  int nconv = geigs.compute(SortRule::LargestMagn);

  // Retrieve results
  Eigen::VectorXd evalues;
  Eigen::MatrixXd evecs;
  if (geigs.info() == CompInfo::Successful) {
    evalues = geigs.eigenvalues();
    evecs = geigs.eigenvectors();
  }

  std::cout << "Number of converged generalized eigenvalues: " << nconv
            << std::endl;
  std::cout << "Generalized eigenvalues found:\n" << evalues << std::endl;
  std::cout << "Generalized eigenvectors found:\n"
            << evecs.topRows(10) << std::endl;
            */
  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////

  return 0;
}