#include <iostream>
#include <PlanetaryModel/All>
#include <new_coupling/Timer>
#include "sem_full.h"
// #include "sem_spheroidal_debug.h"
#include "../SpectraSolver/SpectraSolver/FF"
// #include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/filter_base.h"
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"
#include "read_station.h"
#include "input_parser.h"   // Use the new input parser
#include "read_yspec.h"
#include "read_mineos.h"
#include "SourceInfo.h"
#include "start_element.h"
#include "precon_test.h"
#include <chrono>
#include "spectra_master.h"
#include "ParamInfo.h"
#include <omp.h>
#include "FEM_Preconditioner.h"
#include "BiCGSTABT.h"
#include "ParamRedInfo.h"
#include "SR_Info.h"
// #include <Eigen/UmfPackSupport>

namespace SPARSESPEC {

class Sparse_F_Spec {
public:
  Sparse_F_Spec() {};
  ~Sparse_F_Spec() {};

  auto FrequencySpectrum(SpectraSolver::FreqFull &, Full1D::sem &,
                         SourceInfo::EarthquakeCMT &, InputParameters &,
                         double);

  auto FrequencySpectrum_TEST(SpectraSolver::FreqFull &, Full1D::sem &,
                              SourceInfo::EarthquakeCMT &, InputParameters &,
                              double);

  auto FrequencySpectrum_TEST2(SpectraSolver::FreqFull &, Full1D::sem &,
                               SourceInfo::EarthquakeCMT &, InputParameters &,
                               double);
  template <class model1d>
  auto FrequencySpectrum_TEST_CLEAN(SpectraSolver::FreqFull &, Full1D::sem &,
                                    model1d &, SourceInfo::EarthquakeCMT &,
                                    InputParameters &, double);

  auto FrequencySpectrum_TEST_OUTPUT(SpectraSolver::FreqFull &, Full1D::sem &,
                                     SourceInfo::EarthquakeCMT &,
                                     InputParameters &, double);

  template <class model1d>
  auto FrequencySpectrum_TEST_SPECSEM(SpectraSolver::FreqFull &,
                                      Full1D::specsem &, model1d &,
                                      SourceInfo::EarthquakeCMT &,
                                      InputParameters &, int = 10);

  template <class model1d>
  auto FrequencySpectrum_Variable(SpectraSolver::FreqFull &, model1d &,
                                  SourceInfo::EarthquakeCMT &,
                                  InputParameters &, int, int, int);

  template <class model1d>
  auto FrequencySpectrum_RED(SpectraSolver::FreqFull &, model1d &,
                             SourceInfo::EarthquakeCMT &, InputParameters &,
                             int, int, int, SRInfo &);

private:
};

auto
Sparse_F_Spec::FrequencySpectrum(SpectraSolver::FreqFull &myff,
                                 Full1D::sem &sem,
                                 SourceInfo::EarthquakeCMT &cmt,
                                 InputParameters &params, double droptol) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver;
  // Eigen::BiCGSTAB<SMATRIX, Eigen::IncompleteLUT<Complex>> bicgstab_solver;
  // Eigen::BiCGSTAB<SMATRIX, Eigen::LUPD<Complex>> bicgstab_solver;
  Eigen::BiCGSTAB<SMATRIX, Eigen::DiagonalPreconditioner<Complex>>
      bicgstab_solver;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    std::cout << "Doing Radial Modes\n";
    SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
    SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();

    // calculate force vector
    MATRIX f_r = sem.CalculateForce_R(cmt);

    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr));
    }

    // iterate over frequencies
    for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
      // complex frequency
      Complex w = vec_w[idx] + ieps;

      // force vector at frequency
      MATRIX vec_fw = f_r / (myi * w);

      // build matrix and solve
      SMATRIX w_r = -w * w * in_r + ke_r;
      w_r.makeCompressed();
      solver.compute(w_r);
      MATRIX vec_x = solver.solve(vec_fw);

      // compute responses
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

        // receiver vectors
        // auto RV_Z = sem.RV_Z_R(params, idxr);

        // index
        auto idxpl = 3 * idxr;

        // find response
        vec_raw(idxpl, idx) -= w * w * vec_RV_Z[idxr].cwiseProduct(vec_x).sum();
      }
    };
  }

  //   do toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
      std::cout << "Toroidal,  l = " << idxl << "\n";
      timer1.start();
      SMATRIX mat_ke = sem.MAT_KE_T(idxl).cast<Complex>();
      SMATRIX mat_inertia = sem.MAT_IN_T(idxl).cast<Complex>();
      mat_ke.makeCompressed();
      mat_inertia.makeCompressed();
      timer1.stop("Time for matrix assembling");

      timer1.start();
      // calculate force vector
      MATRIX vec_force = sem.CalculateForce_T(cmt, idxl);
      timer1.stop("Time for force calculation");

      std::vector<MATRIX> vec_RV_THETA, vec_RV_PHI;
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        vec_RV_THETA.push_back(sem.RV_THETA_T(params, idxl, idxr));
        vec_RV_PHI.push_back(sem.RV_PHI_T(params, idxl, idxr));
      }

      timer1.start();
      // iterate over frequencies
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        // complex frequency
        Complex w = vec_w[idx] + ieps;

        // force vector at frequency
        MATRIX vec_fw = vec_force / (myi * w);

        // build matrix and solve
        SMATRIX mat_w = -w * w * mat_inertia + mat_ke;

        // compress, decompose and solve
        mat_w.makeCompressed();
        solver.compute(mat_w);
        MATRIX vec_x = solver.solve(vec_fw);

        // compute responses
        for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

          // receiver vectors
          // auto RV_THETA = sem.RV_THETA_T(params, idxl, idxr);
          // auto RV_PHI = sem.RV_PHI_T(params, idxl, idxr);

          // index
          auto idxpl = 3 * idxr;

          // find response
          vec_raw(idxpl + 1, idx) +=
              w * w * vec_RV_THETA[idxr].cwiseProduct(vec_x).sum();
          vec_raw(idxpl + 2, idx) -=
              w * w * vec_RV_PHI[idxr].cwiseProduct(vec_x).sum();

          // NOTE: signs are according to the convention used in YSpec. In
          // particular for theta its reversed so it is the north componenet
        }
      };
      timer1.stop("Time for frequency loop");
      std::cout << "\n";
    }
  }

  // do spheroidals:
  if (inc_sph) {
    std::cout << "Doing Spheroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1) + 1;
    for (auto idx : rec_elems) {
      std::cout << "Receiver element: " << idx << "\n";
    }
    int lenidx = upidx - lowidx;
    std::cout << "Receiver element indices from " << lowidx << " to "
              << upidx - 1 << "\n";
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
      timer1.start();
      auto total_solver_duration = std::chrono::microseconds::zero();
      auto total_factorize_duration = std::chrono::microseconds::zero();
      auto total_receiver_duration = std::chrono::microseconds::zero();
      auto tot_setup = std::chrono::microseconds::zero();
      auto total_makemat_duration = std::chrono::microseconds::zero();

      double numit = 0;
      bicgstab_solver.setTolerance(1e-6);

      std::cout << "Spheroidal,  l = " << idxl << "\n";
      auto start_mat = std::chrono::high_resolution_clock::now();
      // timer1.start();
      SMATRIX ke_s = sem.MAT_KE(idxl).cast<Complex>();
      SMATRIX in_s = sem.MAT_IN(idxl).cast<Complex>();
      ke_s.makeCompressed();
      in_s.makeCompressed();

      // using tripletlist method
      // auto vec_ke_idx = sem.tripletlist_ke_s_idx(idxl);
      // auto vec_in_idx = sem.tripletlist_in_s_idx(idxl);
      // auto vec_ke_val = sem.tripletlist_ke_s_val(idxl);
      // auto vec_in_val = sem.tripletlist_in_s_val(idxl);
      // std::cout << "Here 1\n";
      // timer1.stop("Time for matrix assembling");
      // analyze pattern for quicker lu decomposition
      // timer1.start();
      {
        SMATRIX mat_test = ke_s + myi * in_s;
        solver.analyzePattern(mat_test);
      }
      // timer1.stop("Time for analyze pattern");
      // using T = Eigen::Triplet<Complex>;
      // std::vector<T> tpl_w;
      // auto numke = vec_ke_idx.size();
      // tpl_w.resize(numke + vec_in_idx.size());
      // for (std::size_t k = 0; k < numke; ++k) {
      //   tpl_w[k] = T(vec_ke_idx[k][0], vec_ke_idx[k][1], vec_ke_val[k]);
      // }
      // timer1.start();
      // calculate force vector
      MATRIX vec_force_sph = sem.CalculateForce(cmt, idxl);
      // timer1.stop("Time for force calculation");
      // timer1.start();
      std::vector<MATRIX> vec_RV_Z, vec_RV_THETA, vec_RV_PHI;
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        vec_RV_Z.push_back(sem.RV_Z(params, idxl, idxr));
        vec_RV_THETA.push_back(sem.RV_THETA(params, idxl, idxr));
        vec_RV_PHI.push_back(sem.RV_PHI(params, idxl, idxr));
      }

      // std::cout << std::setprecision(4)
      //           << vec_RV_Z[0].block(lowidx, 0, lenidx, 2 * idxl + 1) <<
      //           "\n";

      // timer1.stop("Time for receiver vector calculation");
      auto stop_mat = std::chrono::high_resolution_clock::now();
      tot_setup += std::chrono::duration_cast<std::chrono::microseconds>(
          stop_mat - start_mat);
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {

        auto t_ss = std::chrono::high_resolution_clock::now();
        // setup
        Complex w = vec_w[idx] + ieps;
        MATRIX vec_fw_sph = vec_force_sph / (myi * w);
        SMATRIX mat_w_sph = ke_s;
        mat_w_sph -= w * w * in_s;
        mat_w_sph.makeCompressed();

        auto t_sts = std::chrono::high_resolution_clock::now();
        total_makemat_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_sts - t_ss);

        auto t_sf = std::chrono::high_resolution_clock::now();
        solver.factorize(mat_w_sph);
        auto t_stf = std::chrono::high_resolution_clock::now();
        total_factorize_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_stf - t_sf);

        auto start_solver = std::chrono::high_resolution_clock::now();
        MATRIX vec_x_sph = solver.solve(vec_fw_sph);
        auto stop_solver = std::chrono::high_resolution_clock::now();
        total_solver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(stop_solver -
                                                                  start_solver);
        // std::cout << "Here 2\n";
        // bicgstab_solver.preconditioner().setDroptol(droptol);
        // bicgstab_solver.compute(mat_w_sph);
        // if (!((idx - myff.i1()) % 10)) {
        //   bicgstab_solver.preconditioner().addmatrix(mat_w_sph);
        // }
        // MATRIX vec_x_sph = bicgstab_solver.solve(vec_fw_sph);
        // numit += bicgstab_solver.iterations();

        // compute responses
        auto start_receiver = std::chrono::high_resolution_clock::now();
        for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

          // receiver vectors
          //   auto RV_Z = sem.RV_Z(params, idxl, idxr);
          //   auto RV_THETA = sem.RV_THETA(params, idxl, idxr);
          //   auto RV_PHI = sem.RV_PHI(params, idxl, idxr);

          // index
          auto idxpl = 3 * idxr;

          // find response
          vec_raw(idxpl, idx) -= w * w *
                                 vec_RV_Z[idxr]
                                     .block(lowidx, 0, lenidx, 2 * idxl + 1)
                                     .cwiseProduct(vec_x_sph.block(
                                         lowidx, 0, lenidx, 2 * idxl + 1))
                                     .sum();
          vec_raw(idxpl + 1, idx) += w * w *
                                     vec_RV_THETA[idxr]
                                         .block(lowidx, 0, lenidx, 2 * idxl + 1)
                                         .cwiseProduct(vec_x_sph.block(
                                             lowidx, 0, lenidx, 2 * idxl + 1))
                                         .sum();
          vec_raw(idxpl + 2, idx) -= w * w *
                                     vec_RV_PHI[idxr]
                                         .block(lowidx, 0, lenidx, 2 * idxl + 1)
                                         .cwiseProduct(vec_x_sph.block(
                                             lowidx, 0, lenidx, 2 * idxl + 1))
                                         .sum();

          // NOTE: signs are according to the convention used in YSpec. In
          // particular for theta its reversed so it is the north componenet
        }
        auto stop_receiver = std::chrono::high_resolution_clock::now();
        total_receiver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(
                stop_receiver - start_receiver);
      };
      // timer1.stop("Time for frequency loop");
      std::cout << "Average # iterations: "
                << numit / static_cast<double>(myff.i2() - myff.i1()) << "\n";
      std::cout << "Total time for setup: " << tot_setup.count() / 1e6
                << " s\n";
      std::cout << "Total time for making matrices: "
                << total_makemat_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for factorization: "
                << total_factorize_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for matrix solver: "
                << total_solver_duration.count() / 1e6 << " s  \n";
      std::cout << "Total time for receiver vector multiplication: "
                << total_receiver_duration.count() / 1e6 << " s\n";
      timer1.stop("Total time for l = " + std::to_string(idxl));

      std::cout << "\n";
    }
  }

  return vec_raw;
};

auto
Sparse_F_Spec::FrequencySpectrum_TEST(SpectraSolver::FreqFull &myff,
                                      Full1D::sem &sem,
                                      SourceInfo::EarthquakeCMT &cmt,
                                      InputParameters &params, double droptol) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver;
  // Eigen::BiCGSTAB<SMATRIX, Eigen::IncompleteLUT<Complex>> bicgstab_solver;
  // Eigen::BiCGSTAB<SMATRIX, Eigen::LUPD<Complex>> bicgstab_solver;
  Eigen::BiCGSTAB<SMATRIX, Eigen::DiagonalPreconditioner<Complex>>
      bicgstab_solver;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    std::cout << "Doing Radial Modes\n";
    SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
    SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();

    // calculate force vector
    MATRIX f_r = sem.CalculateForce_R(cmt);

    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr));
    }

    // iterate over frequencies
    for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
      // complex frequency
      Complex w = vec_w[idx] + ieps;

      // force vector at frequency
      MATRIX vec_fw = f_r / (myi * w);

      // build matrix and solve
      SMATRIX w_r = -w * w * in_r + ke_r;
      w_r.makeCompressed();
      solver.compute(w_r);
      MATRIX vec_x = solver.solve(vec_fw);

      // compute responses
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

        // receiver vectors
        // auto RV_Z = sem.RV_Z_R(params, idxr);

        // index
        auto idxpl = 3 * idxr;

        // find response
        vec_raw(idxpl, idx) -= w * w * vec_RV_Z[idxr].cwiseProduct(vec_x).sum();
      }
    };
  }

  //   do toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
      std::cout << "Toroidal,  l = " << idxl << "\n";
      timer1.start();
      SMATRIX mat_ke = sem.MAT_KE_T(idxl).cast<Complex>();
      SMATRIX mat_inertia = sem.MAT_IN_T(idxl).cast<Complex>();
      mat_ke.makeCompressed();
      mat_inertia.makeCompressed();
      timer1.stop("Time for matrix assembling");

      timer1.start();
      // calculate force vector
      MATRIX vec_force = sem.CalculateForce_T(cmt, idxl);
      timer1.stop("Time for force calculation");

      std::vector<MATRIX> vec_RV_THETA, vec_RV_PHI;
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        vec_RV_THETA.push_back(sem.RV_THETA_T(params, idxl, idxr));
        vec_RV_PHI.push_back(sem.RV_PHI_T(params, idxl, idxr));
      }

      timer1.start();
      // iterate over frequencies
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        // complex frequency
        Complex w = vec_w[idx] + ieps;

        // force vector at frequency
        MATRIX vec_fw = vec_force / (myi * w);

        // build matrix and solve
        SMATRIX mat_w = -w * w * mat_inertia + mat_ke;

        // compress, decompose and solve
        mat_w.makeCompressed();
        solver.compute(mat_w);
        MATRIX vec_x = solver.solve(vec_fw);

        // compute responses
        for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

          // receiver vectors
          // auto RV_THETA = sem.RV_THETA_T(params, idxl, idxr);
          // auto RV_PHI = sem.RV_PHI_T(params, idxl, idxr);

          // index
          auto idxpl = 3 * idxr;

          // find response
          vec_raw(idxpl + 1, idx) +=
              w * w * vec_RV_THETA[idxr].cwiseProduct(vec_x).sum();
          vec_raw(idxpl + 2, idx) -=
              w * w * vec_RV_PHI[idxr].cwiseProduct(vec_x).sum();

          // NOTE: signs are according to the convention used in YSpec. In
          // particular for theta its reversed so it is the north componenet
        }
      };
      timer1.stop("Time for frequency loop");
      std::cout << "\n";
    }
  }

  // do spheroidals:
  if (inc_sph) {
    std::cout << "Doing Spheroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1) + 1;
    for (auto idx : rec_elems) {
      std::cout << "Receiver element: " << idx << "\n";
    }
    int lenidx = upidx - lowidx;
    std::cout << "Receiver element indices from " << lowidx << " to "
              << upidx - 1 << "\n";
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
      timer1.start();
      auto total_solver_duration = std::chrono::microseconds::zero();
      auto total_factorize_duration = std::chrono::microseconds::zero();
      auto total_receiver_duration = std::chrono::microseconds::zero();
      auto tot_setup = std::chrono::microseconds::zero();
      auto total_makemat_duration = std::chrono::microseconds::zero();

      double numit = 0;
      bicgstab_solver.setTolerance(1e-6);

      std::cout << "Spheroidal,  l = " << idxl << "\n";
      auto start_mat = std::chrono::high_resolution_clock::now();
      // timer1.start();
      SMATRIX ke_s = sem.MAT_KE(idxl).cast<Complex>();
      SMATRIX in_s = sem.MAT_IN(idxl).cast<Complex>();
      ke_s.makeCompressed();
      in_s.makeCompressed();

      // using tripletlist method
      // auto vec_ke_idx = sem.tripletlist_ke_s_idx(idxl);
      // auto vec_in_idx = sem.tripletlist_in_s_idx(idxl);
      // auto vec_ke_val = sem.tripletlist_ke_s_val(idxl);
      // auto vec_in_val = sem.tripletlist_in_s_val(idxl);
      // std::cout << "Here 1\n";
      // timer1.stop("Time for matrix assembling");
      // analyze pattern for quicker lu decomposition
      // timer1.start();
      {
        SMATRIX mat_test = ke_s + myi * in_s;
        solver.analyzePattern(mat_test);
      }
      // timer1.stop("Time for analyze pattern");
      // using T = Eigen::Triplet<Complex>;
      // std::vector<T> tpl_w;
      // auto numke = vec_ke_idx.size();
      // tpl_w.resize(numke + vec_in_idx.size());
      // for (std::size_t k = 0; k < numke; ++k) {
      //   tpl_w[k] = T(vec_ke_idx[k][0], vec_ke_idx[k][1], vec_ke_val[k]);
      // }
      // timer1.start();

      // timer1.stop("Time for force calculation");
      // timer1.start();
      // calculate force vector
      MATRIX vec_force_sph = sem.CalculateForce(cmt, idxl);

      std::vector<MATRIX> vec_RV_Z, vec_RV_THETA, vec_RV_PHI, vec_vals,
          vec_base, vec_forces;
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        vec_RV_Z.push_back(sem.RV_Z(params, idxl, idxr));
        vec_RV_THETA.push_back(sem.RV_THETA(params, idxl, idxr));
        vec_RV_PHI.push_back(sem.RV_PHI(params, idxl, idxr));
        vec_base.push_back(sem.RV_BASE_Z(params, idxl, idxr));

        vec_vals.push_back(sem.RV_VAL_Z(params, idxl, idxr));
        vec_forces.push_back(vec_force_sph * vec_vals.back());
      }

      MATRIX mat_vals(2 * idxl + 1, 3 * params.num_receivers());
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        mat_vals.col(idxr) = sem.RV_VAL_Z(params, idxl, idxr);
        mat_vals.col(idxr + params.num_receivers()) =
            sem.RV_VAL_THETA(params, idxl, idxr);
        mat_vals.col(idxr + 2 * params.num_receivers()) =
            sem.RV_VAL_PHI(params, idxl, idxr);
      }

      MATRIX mat_forces(vec_force_sph.rows(), 3 * params.num_receivers());
      mat_forces = vec_force_sph * mat_vals;   // precompute all force vectors

      MATRIX RV_BASE(2 * params.num_receivers(), vec_force_sph.rows());
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        RV_BASE.row(idxr) = sem.RV_BASE_Z(params, idxl, idxr);
        RV_BASE.row(idxr + params.num_receivers()) =
            sem.RV_BASE_THETA(params, idxl, idxr);
      }

      // MATRIX vec_force_test = vec_force_sph * vec_vals[0];
      // MATRIX vec_force_test2 = vec_force_sph * vec_vals[1];
      // std::cout << std::setprecision(4)
      //           << vec_RV_Z[0].block(lowidx, 0, lenidx, 2 * idxl + 1) <<
      //           "\n";

      // timer1.stop("Time for receiver vector calculation");
      auto stop_mat = std::chrono::high_resolution_clock::now();
      tot_setup += std::chrono::duration_cast<std::chrono::microseconds>(
          stop_mat - start_mat);
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {

        auto t_ss = std::chrono::high_resolution_clock::now();
        // setup
        Complex w = vec_w[idx] + ieps;
        MATRIX vec_fw_sph = vec_force_sph / (myi * w);
        SMATRIX mat_w_sph = ke_s;
        mat_w_sph -= w * w * in_s;
        mat_w_sph.makeCompressed();

        auto t_sts = std::chrono::high_resolution_clock::now();
        total_makemat_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_sts - t_ss);

        auto t_sf = std::chrono::high_resolution_clock::now();
        solver.factorize(mat_w_sph);
        auto t_stf = std::chrono::high_resolution_clock::now();
        total_factorize_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_stf - t_sf);

        auto start_solver = std::chrono::high_resolution_clock::now();
        // MATRIX vec_x_sph = solver.solve(vec_fw_sph);
        MATRIX vec_x_base_sol = solver.solve(mat_forces);
        auto stop_solver = std::chrono::high_resolution_clock::now();
        total_solver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(stop_solver -
                                                                  start_solver);
        // std::cout << "Here 2\n";
        // bicgstab_solver.preconditioner().setDroptol(droptol);
        // bicgstab_solver.compute(mat_w_sph);
        // if (!((idx - myff.i1()) % 10)) {
        //   bicgstab_solver.preconditioner().addmatrix(mat_w_sph);
        // }
        // MATRIX vec_x_sph = bicgstab_solver.solve(vec_fw_sph);
        // numit += bicgstab_solver.iterations();

        // compute responses
        auto start_receiver = std::chrono::high_resolution_clock::now();
        for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

          // receiver vectors
          //   auto RV_Z = sem.RV_Z(params, idxl, idxr);
          //   auto RV_THETA = sem.RV_THETA(params, idxl, idxr);
          //   auto RV_PHI = sem.RV_PHI(params, idxl, idxr);

          // index
          auto idxpl = 3 * idxr;
          auto idxz = idxr;
          auto idxtheta = idxr + params.num_receivers();
          auto idxphi = idxr + 2 * params.num_receivers();
          // find response
          // vec_raw(idxpl, idx) -= w * w *
          //                        vec_RV_Z[idxr]
          //                            .block(lowidx, 0, lenidx, 2 * idxl + 1)
          //                            .cwiseProduct(vec_x_sph.block(
          //                                lowidx, 0, lenidx, 2 * idxl + 1))
          //                            .sum();
          vec_raw(idxpl, idx) -=
              w / myi * (RV_BASE.row(idxz) * vec_x_base_sol.col(idxz)).sum();
          // vec_raw(idxpl + 1, idx) += w * w *
          //                            vec_RV_THETA[idxr]
          //                                .block(lowidx, 0, lenidx, 2 * idxl +
          //                                1) .cwiseProduct(vec_x_sph.block(
          //                                    lowidx, 0, lenidx, 2 * idxl +
          //                                    1))
          //                                .sum();
          vec_raw(idxpl + 1, idx) +=
              w / myi *
              (RV_BASE.row(idxtheta) * vec_x_base_sol.col(idxtheta)).sum();
          vec_raw(idxpl + 2, idx) -=
              w / myi *
              (RV_BASE.row(idxtheta) * vec_x_base_sol.col(idxphi)).sum();
          // vec_raw(idxpl + 2, idx) -= w * w *
          //                            vec_RV_PHI[idxr]
          //                                .block(lowidx, 0, lenidx, 2 * idxl +
          //                                1) .cwiseProduct(vec_x_sph.block(
          //                                    lowidx, 0, lenidx, 2 * idxl +
          //                                    1))
          //                                .sum();

          // NOTE: signs are according to the convention used in YSpec. In
          // particular for theta its reversed so it is the north componenet
        }
        auto stop_receiver = std::chrono::high_resolution_clock::now();
        total_receiver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(
                stop_receiver - start_receiver);
      };
      // timer1.stop("Time for frequency loop");
      std::cout << "Average # iterations: "
                << numit / static_cast<double>(myff.i2() - myff.i1()) << "\n";
      std::cout << "Total time for setup: " << tot_setup.count() / 1e6
                << " s\n";
      std::cout << "Total time for making matrices: "
                << total_makemat_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for factorization: "
                << total_factorize_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for matrix solver: "
                << total_solver_duration.count() / 1e6 << " s  \n";
      std::cout << "Total time for receiver vector multiplication: "
                << total_receiver_duration.count() / 1e6 << " s\n";
      timer1.stop("Total time for l = " + std::to_string(idxl));

      std::cout << "\n";
    }
  }

  return vec_raw;
};

auto
Sparse_F_Spec::FrequencySpectrum_TEST2(SpectraSolver::FreqFull &myff,
                                       Full1D::sem &sem,
                                       SourceInfo::EarthquakeCMT &cmt,
                                       InputParameters &params,
                                       double droptol) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver;
  // Eigen::BiCGSTAB<SMATRIX, Eigen::IncompleteLUT<Complex>> bicgstab_solver;
  // Eigen::BiCGSTAB<SMATRIX, Eigen::LUPD<Complex>> bicgstab_solver;
  Eigen::BiCGSTAB<SMATRIX, Eigen::DiagonalPreconditioner<Complex>>
      bicgstab_solver;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    timer1.start();

    std::cout << "Doing Radial Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
    auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
    for (auto idx : rec_elems) {
      std::cout << "Receiver element: " << idx << "\n";
    }
    int lenidx = upidx - lowidx + 1;
    std::cout << "Receiver element indices from " << lowidx << " to " << upidx
              << "\n";
    SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
    SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();

    // calculate force vector
    MATRIX f_r = sem.CalculateForce_R(cmt);

    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr).block(lowidx, 0, lenidx, 1));
    }

    // iterate over frequencies
    for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
      // complex frequency
      Complex w = vec_w[idx] + ieps;

      // force vector at frequency
      // MATRIX vec_fw = f_r / (myi * w);

      // build matrix and solve
      SMATRIX w_r = -w * w * in_r + ke_r;
      w_r.makeCompressed();
      solver.compute(w_r);
      MATRIX vec_x = solver.solve(f_r);

      // compute responses
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

        // receiver vectors
        // auto RV_Z = sem.RV_Z_R(params, idxr);

        // index
        auto idxpl = 3 * idxr;

        // find response
        vec_raw(idxpl, idx) -=
            w / myi *
            vec_RV_Z[idxr]
                .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                .sum();
      }
    };
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_T(rec_elems[0], 0);
    auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
    for (auto idx : rec_elems) {
      std::cout << "Receiver element: " << idx << "\n";
    }
    int lenidx = upidx - lowidx;
    std::cout << "Receiver element indices from " << lowidx << " to "
              << upidx - 1 << "\n";

    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
      timer1.start();
      auto total_solver_duration = std::chrono::microseconds::zero();
      auto total_factorize_duration = std::chrono::microseconds::zero();
      auto total_receiver_duration = std::chrono::microseconds::zero();
      auto tot_setup = std::chrono::microseconds::zero();
      auto total_makemat_duration = std::chrono::microseconds::zero();

      // double numit = 0;
      // bicgstab_solver.setTolerance(1e-6);

      std::cout << "Toroidal,  l = " << idxl << "\n";
      auto start_mat = std::chrono::high_resolution_clock::now();
      // timer1.start();
      SMATRIX mat_ke_tor = sem.MAT_KE_T(idxl).cast<Complex>();
      SMATRIX mat_in_tor = sem.MAT_IN_T(idxl).cast<Complex>();
      mat_ke_tor.makeCompressed();
      mat_in_tor.makeCompressed();

      {
        SMATRIX mat_test = mat_ke_tor + myi * mat_in_tor;
        solver.analyzePattern(mat_test);
      }

      // calculate force vector
      MATRIX vec_force_tor = sem.CalculateForce_T(cmt, idxl);

      // matrix of receiver vectors
      auto num_r = params.num_receivers();
      MATRIX RV_VALS = MATRIX::Zero(3 * num_r, 2 * idxl + 1);

      for (int idxr = 0; idxr < num_r; ++idxr) {
        RV_VALS.row(idxr * 3 + 1) =
            -sem.RV_VAL_THETA_T(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 2) =
            sem.RV_VAL_PHI_T(params, idxl, idxr).transpose();
      }

      MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);

      MATRIX RED_C = RV_VALS * F_VALS;

      MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);

      MATRIX mat_vals = MATRIX::Zero(2 * idxl + 1, 3 * num_r);
      for (int idxr = 0; idxr < num_r; ++idxr) {
        // mat_vals.col(idxr * 3) = sem.RV_VAL_Z_T(params, idxl, idxr);
        mat_vals.col(idxr * 3 + 1) = sem.RV_VAL_THETA_T(params, idxl, idxr);
        mat_vals.col(idxr * 3 + 2) = sem.RV_VAL_PHI_T(params, idxl, idxr);
      };

      MATRIX mat_forces(vec_force_tor.rows(), 3 * num_r);
      mat_forces = vec_force_tor * mat_vals;   // precompute all force vectors

      auto nrec = num_r;

      MATRIX RV_BASE = MATRIX::Zero(3 * nrec, lenidx);
      for (int idxr = 0; idxr < nrec; ++idxr) {
        RV_BASE.row(idxr * 3 + 1) =
            sem.RV_BASE_THETA_T(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 2) =
            sem.RV_BASE_PHI_T(params, idxl, idxr).block(0, lowidx, 1, lenidx);
      }

      auto stop_mat = std::chrono::high_resolution_clock::now();
      tot_setup += std::chrono::duration_cast<std::chrono::microseconds>(
          stop_mat - start_mat);
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {

        ////////////////////////////////////////////////////////////////////////
        auto t_ss = std::chrono::high_resolution_clock::now();

        // setup
        Complex w = vec_w[idx] + ieps;
        SMATRIX mat_w_tor = mat_ke_tor - w * w * mat_in_tor;
        mat_w_tor.makeCompressed();

        auto t_sts = std::chrono::high_resolution_clock::now();
        total_makemat_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_sts - t_ss);

        ////////////////////////////////////////////////////////////////////////

        auto t_sf = std::chrono::high_resolution_clock::now();
        solver.factorize(mat_w_tor);
        auto t_stf = std::chrono::high_resolution_clock::now();
        total_factorize_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_stf - t_sf);

        ////////////////////////////////////////////////////////////////////////
        auto start_solver = std::chrono::high_resolution_clock::now();

        MATRIX vec_coeff_sol = solver.solve(F_BASE);
        auto stop_solver = std::chrono::high_resolution_clock::now();
        total_solver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(stop_solver -
                                                                  start_solver);
        ////////////////////////////////////////////////////////////////////////

        // compute responses
        auto start_receiver = std::chrono::high_resolution_clock::now();
        // std::cout << "Size: " << vec_coeff_sol.rows() << " "
        //           << vec_coeff_sol.cols() << "\n";
        auto testmult = RV_BASE * vec_coeff_sol.block(lowidx, 0, lenidx, 2);
        vec_raw.col(idx) -=
            w / myi * RED_C.cwiseProduct(testmult).rowwise().sum();

        auto stop_receiver = std::chrono::high_resolution_clock::now();
        total_receiver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(
                stop_receiver - start_receiver);
      };
      timer1.stop("Time for frequency loop");
      // std::cout << "Average # iterations: "
      //           << numit / static_cast<double>(myff.i2() - myff.i1()) <<
      //           "\n";
      std::cout << "Total time for setup: " << tot_setup.count() / 1e6
                << " s\n";
      std::cout << "Total time for making matrices: "
                << total_makemat_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for factorization: "
                << total_factorize_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for matrix solver: "
                << total_solver_duration.count() / 1e6 << " s  \n";
      std::cout << "Total time for receiver vector multiplication: "
                << total_receiver_duration.count() / 1e6 << " s\n";
      timer1.stop("Total time for l = " + std::to_string(idxl));

      // std::cout << "\n";
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // do spheroidals:
  if (inc_sph) {
    std::cout << "Doing Spheroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
    for (auto idx : rec_elems) {
      std::cout << "Receiver element: " << idx << "\n";
    }
    int lenidx = upidx - lowidx + 1;
    std::cout << "Receiver element indices from " << lowidx << " to " << upidx
              << "\n";
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
      timer1.start();
      auto total_solver_duration = std::chrono::microseconds::zero();
      auto total_factorize_duration = std::chrono::microseconds::zero();
      auto total_receiver_duration = std::chrono::microseconds::zero();
      auto tot_setup = std::chrono::microseconds::zero();
      auto total_makemat_duration = std::chrono::microseconds::zero();

      // double numit = 0;
      // bicgstab_solver.setTolerance(1e-6);

      std::cout << "Spheroidal,  l = " << idxl << "\n";
      auto start_mat = std::chrono::high_resolution_clock::now();
      // timer1.start();
      SMATRIX ke_s = sem.MAT_KE(idxl).cast<Complex>();
      SMATRIX in_s = sem.MAT_IN(idxl).cast<Complex>();
      ke_s.makeCompressed();
      in_s.makeCompressed();

      {
        SMATRIX mat_test = ke_s + myi * in_s;
        solver.analyzePattern(mat_test);
      }

      // calculate force vector
      MATRIX vec_force_sph = sem.CalculateForce(cmt, idxl);

      // matrix of receiver vectors
      MATRIX RV_VALS(3 * params.num_receivers(), 2 * idxl + 1);
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        RV_VALS.row(idxr * 3) = sem.RV_VAL_Z(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 1) =
            -sem.RV_VAL_THETA(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 2) =
            sem.RV_VAL_PHI(params, idxl, idxr).transpose();
      }

      MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);

      MATRIX RED_C = RV_VALS * F_VALS;

      MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);

      MATRIX mat_vals(2 * idxl + 1, 3 * params.num_receivers());
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
        mat_vals.col(idxr * 3) = sem.RV_VAL_Z(params, idxl, idxr);
        mat_vals.col(idxr * 3 + 1) = sem.RV_VAL_THETA(params, idxl, idxr);
        mat_vals.col(idxr * 3 + 2) = sem.RV_VAL_PHI(params, idxl, idxr);
      };

      MATRIX mat_forces(vec_force_sph.rows(), 3 * params.num_receivers());
      mat_forces = vec_force_sph * mat_vals;   // precompute all force vectors

      auto nrec = params.num_receivers();

      MATRIX RV_BASE(3 * nrec, lenidx);
      // std::cout << "Debug test 0\n";
      for (int idxr = 0; idxr < nrec; ++idxr) {
        RV_BASE.row(idxr * 3) =
            sem.RV_BASE_Z(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 1) =
            sem.RV_BASE_THETA(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 2) =
            sem.RV_BASE_THETA(params, idxl, idxr).block(0, lowidx, 1, lenidx);
      }
      // std::cout << "Debug test 1\n";

      auto stop_mat = std::chrono::high_resolution_clock::now();
      tot_setup += std::chrono::duration_cast<std::chrono::microseconds>(
          stop_mat - start_mat);
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {

        ////////////////////////////////////////////////////////////////////////
        auto t_ss = std::chrono::high_resolution_clock::now();

        // setup
        Complex w = vec_w[idx] + ieps;
        SMATRIX mat_w_sph = ke_s - w * w * in_s;
        mat_w_sph.makeCompressed();

        auto t_sts = std::chrono::high_resolution_clock::now();
        total_makemat_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_sts - t_ss);

        ////////////////////////////////////////////////////////////////////////

        auto t_sf = std::chrono::high_resolution_clock::now();
        solver.factorize(mat_w_sph);
        auto t_stf = std::chrono::high_resolution_clock::now();
        total_factorize_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_stf - t_sf);

        ////////////////////////////////////////////////////////////////////////
        auto start_solver = std::chrono::high_resolution_clock::now();

        MATRIX vec_coeff_sol = solver.solve(F_BASE);
        auto stop_solver = std::chrono::high_resolution_clock::now();
        total_solver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(stop_solver -
                                                                  start_solver);
        ////////////////////////////////////////////////////////////////////////

        // compute responses
        auto start_receiver = std::chrono::high_resolution_clock::now();
        // std::cout << "Debug test 2\n";
        // std::cout << "lenidx: " << lenidx << ", lowidx: " << lowidx
        //           << ", upidx: " << upidx << "\n";
        // std::cout << "vec_coeff_sol size: " << vec_coeff_sol.rows() << " x "
        //           << vec_coeff_sol.cols() << "\n";
        auto testmult = RV_BASE * vec_coeff_sol.block(lowidx, 0, lenidx, 4);
        // std::cout << "Debug test 3\n";
        vec_raw.col(idx) -=
            w / myi * RED_C.cwiseProduct(testmult).rowwise().sum();

        auto stop_receiver = std::chrono::high_resolution_clock::now();
        total_receiver_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(
                stop_receiver - start_receiver);
      };
      timer1.stop("Time for frequency loop");
      // std::cout << "Average # iterations: "
      //           << numit / static_cast<double>(myff.i2() - myff.i1()) <<
      //           "\n";
      std::cout << "Total time for setup: " << tot_setup.count() / 1e6
                << " s\n";
      std::cout << "Total time for making matrices: "
                << total_makemat_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for factorization: "
                << total_factorize_duration.count() / 1e6 << " s\n";
      std::cout << "Total time for matrix solver: "
                << total_solver_duration.count() / 1e6 << " s  \n";
      std::cout << "Total time for receiver vector multiplication: "
                << total_receiver_duration.count() / 1e6 << " s\n";
      timer1.stop("Total time for l = " + std::to_string(idxl));

      std::cout << "\n";
    }
  }

  return vec_raw;
};

template <class model1d>
auto
Sparse_F_Spec::FrequencySpectrum_TEST_CLEAN(
    SpectraSolver::FreqFull &myff, Full1D::sem &sem, model1d &inp_model,
    SourceInfo::EarthquakeCMT &cmt, InputParameters &params, double droptol) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver, solver1;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();
  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    timer1.start();

    std::cout << "Doing Radial Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
    auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
    SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();

    // calculate force vector
    MATRIX f_r = sem.CalculateForce_R(cmt);

    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < num_rec; ++idxr) {
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr).block(lowidx, 0, lenidx, 1));
    }

    // iterate over frequencies
    for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
      // complex frequency
      Complex w = vec_w[idx] + ieps;

      // build matrix and solve
      SMATRIX w_r = -w * w * in_r + ke_r;
      w_r.makeCompressed();
      solver.compute(w_r);
      MATRIX vec_x = solver.solve(f_r);

      auto twopi = 2.0 * 3.14159265358979323846;

      // compute responses
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        // index
        auto idxpl = 3 * idxr;

        // find response
        vec_raw(idxpl, idx) -=
            w / myi *
            vec_RV_Z[idxr]
                .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                .sum();
      }
    };
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_T(rec_elems[0], 0);
    auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;
    auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;

    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

      // auto total_receiver_duration = std::chrono::microseconds::zero();
      // auto tot_setup = std::chrono::microseconds::zero();
      // auto total_makemat_duration = std::chrono::microseconds::zero();

      std::cout << "Toroidal,  l = " << idxl << "\n";
      SMATRIX mat_ke_tor = sem.MAT_KE_T(idxl).cast<Complex>();
      SMATRIX mat_in_tor = sem.MAT_IN_T(idxl).cast<Complex>();
      mat_ke_tor.makeCompressed();
      mat_in_tor.makeCompressed();
      // analyze pattern
      {
        SMATRIX mat_test = mat_ke_tor + myi * mat_in_tor;
        solver.analyzePattern(mat_test);
      }

      // matrix of receiver vectors
      MATRIX RV_VALS = MATRIX::Zero(3 * num_rec, 2 * idxl + 1);

      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_VALS.row(idxr * 3 + 1) =
            -sem.RV_VAL_THETA_T(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 2) =
            sem.RV_VAL_PHI_T(params, idxl, idxr).transpose();
      }

      MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);
      MATRIX RED_C = RV_VALS * F_VALS;
      MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);

      // auto nrec = num_r;

      MATRIX RV_BASE = MATRIX::Zero(3 * num_rec, lenidx);
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_BASE.row(idxr * 3 + 1) =
            sem.RV_BASE_THETA_T(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 2) =
            sem.RV_BASE_PHI_T(params, idxl, idxr).block(0, lowidx, 1, lenidx);
      }

      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {

        // lower element at this frequency
        auto idxlow_e =
            SpectralTools::StartElement_Tor(sem, inp_model, idxl, vec_w[idx]);
        std::size_t ridx = sem.LtG_T(idxlow_e, 0);
        std::size_t len_ms = lentor - ridx;

        // setup
        Complex w = vec_w[idx] + ieps;
        // SMATRIX mat_w_tor = mat_ke_tor - w * w * mat_in_tor;
        // mat_w_tor.makeCompressed();
        // solver.factorize(mat_w_tor);
        // MATRIX vec_coeff_sol = solver.solve(F_BASE);

        // compute responses
        // auto testmult = RV_BASE * vec_coeff_sol.block(lowidx, 0, lenidx, 2);
        // vec_raw.col(idx) -=
        //     w / myi * RED_C.cwiseProduct(testmult).rowwise().sum();
        // std::cout << "Reduced index: " << ridx
        //           << ", reduced matrix length: " << len_ms << "\n";
        // testing reduced
        SMATRIX mat_w_tor_red =
            mat_ke_tor.block(ridx, ridx, len_ms, len_ms) -
            w * w * mat_in_tor.block(ridx, ridx, len_ms, len_ms);
        mat_w_tor_red.makeCompressed();
        auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
        solver1.compute(mat_w_tor_red);
        MATRIX vec_sol = solver1.solve(f_red);

        auto lidx = lowidx - ridx;
        // compute responses
        auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 2);
        vec_raw.col(idx) -=
            w / myi * RED_C.cwiseProduct(testmult_red).rowwise().sum();
      };
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // do spheroidals:
  if (inc_sph) {
    std::cout << "\nDoing Spheroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;

    // loop over l values
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

      auto tot_setup = std::chrono::microseconds::zero();
      auto tot_mat = std::chrono::microseconds::zero();
      auto total_factorize_duration = std::chrono::microseconds::zero();

      timer1.start();
      std::cout << "Spheroidal,  l = " << idxl << "\n";
      SMATRIX ke_s = sem.MAT_KE(idxl).cast<Complex>();
      SMATRIX in_s = sem.MAT_IN(idxl).cast<Complex>();
      ke_s.makeCompressed();
      in_s.makeCompressed();

      // analyze pattern
      {
        SMATRIX mat_test = ke_s + myi * in_s;
        solver.analyzePattern(mat_test);
      }

      // matrix of receiver vectors
      MATRIX RV_VALS(3 * num_rec, 2 * idxl + 1);
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_VALS.row(idxr * 3) = sem.RV_VAL_Z(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 1) =
            -sem.RV_VAL_THETA(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 2) =
            sem.RV_VAL_PHI(params, idxl, idxr).transpose();
      }

      MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);
      MATRIX RED_C = RV_VALS * F_VALS;
      MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);

      MATRIX RV_BASE(3 * num_rec, lenidx);
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_BASE.row(idxr * 3) =
            sem.RV_BASE_Z(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 1) =
            sem.RV_BASE_THETA(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 2) =
            sem.RV_BASE_THETA(params, idxl, idxr).block(0, lowidx, 1, lenidx);
      }

      auto t_ss = std::chrono::high_resolution_clock::now();
      std::vector<int> vec_ridx(myff.i2() - myff.i1(), 0);
      for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
        auto idxn = myff.i2() - idx - 1;
        int idxlow_e;
        if ((idxn % 10) == 0) {
          idxlow_e =
              SpectralTools::StartElement_Sph(sem, inp_model, idxl, vec_w[idx]);
          std::size_t idx_rs = sem.LtG_S(0, idxlow_e, 0);
          vec_ridx[idx - myff.i1()] = idx_rs;
        } else {
          vec_ridx[idx - myff.i1()] = vec_ridx[idx - myff.i1() + 1];
        }
      }
      auto t_sts = std::chrono::high_resolution_clock::now();
      tot_setup +=
          std::chrono::duration_cast<std::chrono::microseconds>(t_sts - t_ss);

      for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {

        t_ss = std::chrono::high_resolution_clock::now();
        std::size_t idx_rs = vec_ridx[idx - myff.i1()];
        std::size_t len_ms = lensph - idx_rs;
        t_sts = std::chrono::high_resolution_clock::now();

        tot_setup +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_sts - t_ss);
        // setup
        Complex w = vec_w[idx] + ieps;
        // SMATRIX mat_w_sph = ke_s - w * w * in_s;
        // mat_w_sph.makeCompressed();
        // solver.factorize(mat_w_sph);
        // MATRIX vec_coeff_sol = solver.solve(F_BASE);

        // compute responses
        // auto testmult = RV_BASE * vec_coeff_sol.block(lowidx, 0, lenidx, 4);
        // vec_raw.col(idx) -=
        //     w / myi * RED_C.cwiseProduct(testmult).rowwise().sum();

        // testing reduced
        SMATRIX mat_sph = ke_s.block(idx_rs, idx_rs, len_ms, len_ms) -
                          w * w * in_s.block(idx_rs, idx_rs, len_ms, len_ms);
        mat_sph.makeCompressed();
        auto f_red = F_BASE.block(idx_rs, 0, len_ms, F_BASE.cols());
        auto t_stm = std::chrono::high_resolution_clock::now();
        tot_mat += std::chrono::duration_cast<std::chrono::microseconds>(t_stm -
                                                                         t_sts);

        auto t_sf = std::chrono::high_resolution_clock::now();
        auto idxn = myff.i2() - idx - 1;
        if ((idxn % 10) == 0) {
          solver1.compute(mat_sph);
        } else {
          solver1.factorize(mat_sph);
        }

        MATRIX vec_sol = solver1.solve(f_red);

        auto lidx = lowidx - idx_rs;
        // compute responses
        auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);
        vec_raw.col(idx) -=
            w / myi * RED_C.cwiseProduct(testmult_red).rowwise().sum();
        auto t_stf = std::chrono::high_resolution_clock::now();
        total_factorize_duration +=
            std::chrono::duration_cast<std::chrono::microseconds>(t_stf - t_sf);
      };
      std::cout << "Total time for setup: " << tot_setup.count() / 1e6
                << " s\n";
      std::cout << "Total time for making matrices: " << tot_mat.count() / 1e6
                << " s\n";
      std::cout << "Total time for factorization: "
                << total_factorize_duration.count() / 1e6 << " s\n";
      timer1.stop("Time for l = " + std::to_string(idxl));
      std::cout << "\n";
    }
  }

  return vec_raw;
};

auto
Sparse_F_Spec::FrequencySpectrum_TEST_OUTPUT(SpectraSolver::FreqFull &myff,
                                             Full1D::sem &sem,
                                             SourceInfo::EarthquakeCMT &cmt,
                                             InputParameters &params,
                                             double droptol) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();
  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    timer1.start();

    std::cout << "Doing Radial Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
    auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
    SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();

    // calculate force vector
    MATRIX f_r = sem.CalculateForce_R(cmt);

    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < num_rec; ++idxr) {
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr).block(lowidx, 0, lenidx, 1));
    }

    // iterate over frequencies
    for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
      // complex frequency
      Complex w = vec_w[idx] + ieps;

      // build matrix and solve
      SMATRIX w_r = -w * w * in_r + ke_r;
      w_r.makeCompressed();
      solver.compute(w_r);
      MATRIX vec_x = solver.solve(f_r);

      auto twopi = 2.0 * 3.14159265358979323846;

      if ((vec_w[idx] * 1000.0 / (930.0676 * twopi) > 24.428) &&
          (vec_w[idx] * 1000.0 / (930.0676 * twopi) < 24.430)) {
        // std::cout << "Size of vec_coeff_sol: " << vec_coeff_sol.rows()
        //           << " x " << vec_coeff_sol.cols() << "\n";
        std::string pathtofile2 = "./work/spheroidal/full_solution_w.out";
        std::ofstream file2(pathtofile2);
        auto vec_coeff_sol = vec_x;
        // Write header (optional)
        // file <<
        // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
        file2.setf(std::ios::fixed);
        file2 << std::setprecision(16);
        // auto nval = 1.0 / prem.TimeNorm();
        for (int idxe = 0; idxe < sem.mesh().NE(); ++idxe) {
          int maxq = NQ;
          for (int idxq = 0; idxq < maxq; ++idxq) {
            auto uidx = sem.LtG_R(0, idxe, idxq);
            // auto vidx = sem.LtG_S(1, idxe, idxq);
            auto widx = sem.LtG_R(1, idxe, idxq);
            file2 << std::setprecision(22) << sem.mesh().NodeRadius(idxe, idxq)
                  << ';' << vec_coeff_sol(uidx, 0).real() << ';'
                  << vec_coeff_sol(uidx, 0).imag() << ';'
                  << std::abs(vec_coeff_sol(uidx, 0)) << ';'
                  << vec_coeff_sol(widx, 0).real() << ';'
                  << vec_coeff_sol(widx, 0).imag() << ';'
                  << std::abs(vec_coeff_sol(widx, 0)) << ';'
                  << vec_coeff_sol(uidx, 1).real() << ';'
                  << vec_coeff_sol(uidx, 1).imag() << ';'
                  << std::abs(vec_coeff_sol(uidx, 1)) << ';'
                  << vec_coeff_sol(widx, 1).real() << ';'
                  << vec_coeff_sol(widx, 1).imag() << ";"
                  << std::abs(vec_coeff_sol(widx, 1)) << '\n';
          }
        }

        file2.close();
        // std::cout << "Debug freq: " << vec_w[idx] * 1000.0/930 << " Hz\n";
      }

      // compute responses
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        // index
        auto idxpl = 3 * idxr;

        // find response
        vec_raw(idxpl, idx) -=
            w / myi *
            vec_RV_Z[idxr]
                .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                .sum();
      }
    };
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_T(rec_elems[0], 0);
    auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;

    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

      std::cout << "Toroidal,  l = " << idxl << "\n";
      SMATRIX mat_ke_tor = sem.MAT_KE_T(idxl).cast<Complex>();
      SMATRIX mat_in_tor = sem.MAT_IN_T(idxl).cast<Complex>();
      mat_ke_tor.makeCompressed();
      mat_in_tor.makeCompressed();
      // analyze pattern
      {
        SMATRIX mat_test = mat_ke_tor + myi * mat_in_tor;
        solver.analyzePattern(mat_test);
      }

      // matrix of receiver vectors
      MATRIX RV_VALS = MATRIX::Zero(3 * num_rec, 2 * idxl + 1);

      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_VALS.row(idxr * 3 + 1) =
            -sem.RV_VAL_THETA_T(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 2) =
            sem.RV_VAL_PHI_T(params, idxl, idxr).transpose();
      }

      MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);
      MATRIX RED_C = RV_VALS * F_VALS;
      MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);

      // auto nrec = num_r;

      MATRIX RV_BASE = MATRIX::Zero(3 * num_rec, lenidx);
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_BASE.row(idxr * 3 + 1) =
            sem.RV_BASE_THETA_T(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 2) =
            sem.RV_BASE_PHI_T(params, idxl, idxr).block(0, lowidx, 1, lenidx);
      }

      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        // setup
        Complex w = vec_w[idx] + ieps;
        SMATRIX mat_w_tor = mat_ke_tor - w * w * mat_in_tor;
        mat_w_tor.makeCompressed();
        solver.factorize(mat_w_tor);
        MATRIX vec_coeff_sol = solver.solve(F_BASE);

        // compute responses
        auto testmult = RV_BASE * vec_coeff_sol.block(lowidx, 0, lenidx, 2);
        vec_raw.col(idx) -=
            w / myi * RED_C.cwiseProduct(testmult).rowwise().sum();
      };
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // do spheroidals:
  if (inc_sph) {
    std::cout << "\nDoing Spheroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;

    // loop over l values
    for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

      std::cout << "Spheroidal,  l = " << idxl << "\n";
      SMATRIX ke_s = sem.MAT_KE(idxl).cast<Complex>();
      SMATRIX in_s = sem.MAT_IN(idxl).cast<Complex>();
      ke_s.makeCompressed();
      in_s.makeCompressed();

      // analyze pattern
      {
        SMATRIX mat_test = ke_s + myi * in_s;
        solver.analyzePattern(mat_test);
      }

      // matrix of receiver vectors
      MATRIX RV_VALS = MATRIX::Zero(3 * num_rec, 2 * idxl + 1);
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_VALS.row(idxr * 3) = sem.RV_VAL_Z(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 1) =
            -sem.RV_VAL_THETA(params, idxl, idxr).transpose();
        RV_VALS.row(idxr * 3 + 2) =
            sem.RV_VAL_PHI(params, idxl, idxr).transpose();
      }

      // std::cout << "Z force coeff: \n"
      //           << sem.RV_VAL_Z(params, idxl, 0) << "\n\n";

      MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);

      // std::cout << "Force coeff matrix: \n" << F_VALS << "\n\n";
      MATRIX RED_C = RV_VALS * F_VALS;

      // std::cout << "Coefficient multiplication matrix: \n"
      //           << RED_C << "\n\n";

      MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);

      MATRIX RV_BASE = MATRIX::Zero(3 * num_rec, lenidx);
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        RV_BASE.row(idxr * 3) =
            sem.RV_BASE_Z(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 1) =
            sem.RV_BASE_THETA(params, idxl, idxr).block(0, lowidx, 1, lenidx);
        RV_BASE.row(idxr * 3 + 2) =
            sem.RV_BASE_THETA(params, idxl, idxr).block(0, lowidx, 1, lenidx);
      }
      // std::cout << "Base matrix: \n" << RV_BASE << "\n\n";

      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {

        // if (idx == (myff.i2() - 1)) {
        //   std::cout << "Debug: last frequency index " << idx << " "
        //             << vec_w[idx] * 1000.0 / (930.0676 * twopi) << "\n";
        // }

        // setup
        Complex w = vec_w[idx] + ieps;
        SMATRIX mat_w_sph = ke_s - w * w * in_s;
        mat_w_sph.makeCompressed();
        solver.factorize(mat_w_sph);
        MATRIX vec_coeff_sol = solver.solve(F_BASE);
        auto twopi = 2.0 * 3.14159265358979323846;
        // compute responses
        auto testmult = RV_BASE * vec_coeff_sol.block(lowidx, 0, lenidx, 4);
        vec_raw.col(idx) -=
            w / myi * RED_C.cwiseProduct(testmult).rowwise().sum();

        /*
        // if ((vec_w[idx] * 1000.0 / (930.0676 * twopi) > 24.428) &&
        //     (vec_w[idx] * 1000.0 / (930.0676 * twopi) < 24.430)) {

        if ((vec_w[idx] * 1000.0 / (930.0676 * twopi) > 10.397) &&
            (vec_w[idx] * 1000.0 / (930.0676 * twopi) < 10.400)) {
          // std::cout << "Size of vec_coeff_sol: " << vec_coeff_sol.rows()
          //           << " x " << vec_coeff_sol.cols() << "\n";
          std::string pathtofile2 = "./work/spheroidal/full_solution_w.out";
          std::ofstream file2(pathtofile2);

          // Write header (optional)
          // file <<
          //
        "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
          file2.setf(std::ios::fixed);
          file2 << std::setprecision(16);
          // auto nval = 1.0 / prem.TimeNorm();
          for (int idxe = 0; idxe < sem.mesh().NE(); ++idxe) {
            int maxq = NQ;
            // if (idxe == (sem.mesh().NE() - 1)) {
            //   maxq = NQ;
            // }
            // std::cout << "Element: " << idxe
            //           << ", lower radius: " << sem.mesh().NodeRadius(idxe, 0)
            //           << ", upper radius: "
            //           << sem.mesh().NodeRadius(idxe, maxq - 1) << "\n";
            for (int idxq = 0; idxq < maxq; ++idxq) {
              auto uidx = sem.LtG_S(0, idxe, idxq);
              auto vidx = sem.LtG_S(1, idxe, idxq);
              auto widx = sem.LtG_S(2, idxe, idxq);
              file2 << std::setprecision(22)
                    << sem.mesh().NodeRadius(idxe, idxq) << ';'
                    << vec_coeff_sol(uidx, 0).real() << ';'
                    << vec_coeff_sol(uidx, 0).imag() << ';'
                    << std::abs(vec_coeff_sol(uidx, 0)) << ';'
                    << vec_coeff_sol(vidx, 0).real() << ';'
                    << vec_coeff_sol(vidx, 0).imag() << ';'
                    << std::abs(vec_coeff_sol(vidx, 0)) << ';'
                    << vec_coeff_sol(widx, 0).real() << ';'
                    << vec_coeff_sol(widx, 0).imag() << ';'
                    << std::abs(vec_coeff_sol(widx, 0)) << ';'
                    << vec_coeff_sol(uidx, 1).real() << ';'
                    << vec_coeff_sol(uidx, 1).imag() << ';'
                    << std::abs(vec_coeff_sol(uidx, 1)) << ';'
                    << vec_coeff_sol(vidx, 1).real() << ';'
                    << vec_coeff_sol(vidx, 1).imag() << ';'
                    << std::abs(vec_coeff_sol(vidx, 1)) << ';'
                    << vec_coeff_sol(widx, 1).real() << ';'
                    << vec_coeff_sol(widx, 1).imag() << ';'
                    << std::abs(vec_coeff_sol(widx, 1)) << ';'
                    << vec_coeff_sol(uidx, 2).real() << ';'
                    << vec_coeff_sol(uidx, 2).imag() << ';'
                    << std::abs(vec_coeff_sol(uidx, 2)) << ';'
                    << vec_coeff_sol(vidx, 2).real() << ';'
                    << vec_coeff_sol(vidx, 2).imag() << ';'
                    << std::abs(vec_coeff_sol(vidx, 2)) << ';'
                    << vec_coeff_sol(widx, 2).real() << ';'
                    << vec_coeff_sol(widx, 2).imag() << ';'
                    << std::abs(vec_coeff_sol(widx, 2)) << ';'
                    << vec_coeff_sol(uidx, 3).real() << ';'
                    << vec_coeff_sol(uidx, 3).imag() << ';'
                    << std::abs(vec_coeff_sol(uidx, 3)) << ';'
                    << vec_coeff_sol(vidx, 3).real() << ';'
                    << vec_coeff_sol(vidx, 3).imag() << ';'
                    << std::abs(vec_coeff_sol(vidx, 3)) << ';'
                    << vec_coeff_sol(widx, 3).real() << ';'
                    << vec_coeff_sol(widx, 3).imag() << ';'
                    << std::abs(vec_coeff_sol(widx, 3)) << '\n';
            }
          }

          file2.close();
          // std::cout << "Debug freq: " << vec_w[idx] * 1000.0/930 << " Hz\n";

        }*/
      };
    }
  }

  return vec_raw;
};

template <class model1d>
auto
Sparse_F_Spec::FrequencySpectrum_TEST_SPECSEM(
    SpectraSolver::FreqFull &myff, Full1D::specsem &sem, model1d &inp_model,
    SourceInfo::EarthquakeCMT &cmt, InputParameters &params, int nskip) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;
  std::cout << "Running FrequencySpectrum_TEST_SPECSEM with nskip = " << nskip
            << "\n";

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  auto tref = inp_model.TREF();
  auto w0 = 2.0 * 3.14159265358979323846 / tref;
  auto twodivpi = 2.0 / 3.14159265358979323846;

  // std::cout << "Reference frequency: " << w0 << " Hz\n";
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver, solver1;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamInfo param_info(params, lmax);
  auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();

  // std::cout << "Output type: " << params.output_type() << "\n";
  ///////////////////////////////////
  auto idx_source = sem.Source_Element(cmt);
  //   timer1.start();

  // do radials
  if (inc_rad) {
    timer1.start();

    std::cout << "Doing Radial Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
    auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
    SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
    SMATRIX ke_r_atten = sem.MAT_KE_R_ATTEN().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();
    ke_r_atten.makeCompressed();

    auto twopi = 2.0 * 3.14159265358979323846;

    // calculate force vector
    MATRIX f_r = sem.CalculateForce_R(cmt);

    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < num_rec; ++idxr) {
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr).block(lowidx, 0, lenidx, 1));
    }

#pragma omp parallel default(shared) private(solver)
    {
#pragma omp for schedule(dynamic, 10)
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        // complex frequency
        Complex w = vec_w[idx] + ieps;

        // build matrix and solve
        SMATRIX w_r = -w * w * in_r + ke_r;
        // if including attenuation
        if (params.attenuation()) {
          // std::cout << "Including attenuation in toroidal modes\n";
          w_r += (myi + twodivpi * std::log(vec_w[idx] / w0)) * ke_r_atten;
        }
        w_r.makeCompressed();
        solver.compute(w_r);
        MATRIX vec_x = solver.solve(f_r);

        // compute responses
        for (int idxr = 0; idxr < num_rec; ++idxr) {
          // index
          auto idxpl = 3 * idxr;
          Complex mfact;
          if (params.output_type() == 0) {
            mfact = -myi / w;
          } else if (params.output_type() == 1) {
            mfact = 1.0;
          } else if (params.output_type() == 2) {
            mfact = myi * w;
          }
// find response
#pragma omp critical(torvecadd)
          {
            vec_raw(idxpl, idx) +=
                mfact * vec_RV_Z[idxr]
                            .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                            .sum();
          }
        }
      };
    }
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_T(rec_elems[0], 0);
    auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;
    auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

        // auto total_receiver_duration = std::chrono::microseconds::zero();
        // auto tot_setup = std::chrono::microseconds::zero();
        // auto total_makemat_duration = std::chrono::microseconds::zero();

#pragma omp critical(toroutput)
        {
          // std::cout << "Toroidal,  l = " << idxl << "\n";
        }

        // get matrices
        SMATRIX mat_ke_tor = sem.MAT_KE_T_K(idxl).cast<Complex>();
        SMATRIX mat_in_tor = sem.MAT_IN_T_K(idxl).cast<Complex>();
        SMATRIX mat_ke_tor_atten = sem.MAT_KE_T_ATTEN(idxl).cast<Complex>();
        mat_ke_tor.makeCompressed();
        mat_in_tor.makeCompressed();
        mat_ke_tor_atten.makeCompressed();

        // std::cout << "Test 1\n";

        ///////////////////////////////////////////////////////////////////////////
        // setting up forces and receiver vectors
        // MATRIX RV_VALS = sem.RV_FULL_T(params, idxl);
        MATRIX RV_VALS = param_info.RV_FULL_TOR(idxl);
        // std::cout << "Test 2\n";
        MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);
        // std::cout << "Test 3\n";
        MATRIX RED_C = RV_VALS * F_VALS;
        // std::cout << "Test 4\n";
        MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
        // std::cout << "Test 5\n";
        MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);
        // std::cout << "Test 6\n";

        // indices of reduced matrices at each frequency
        auto ridxtor = SpectralTools::AllIndices_TOR(sem, inp_model, idxl, myff,
                                                     idx_source, nskip);
        // std::cout << "Test 7\n";
        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {

          // lower element at this frequency
          // auto idxlow_e =
          //     SpectralTools::StartElement_Tor(sem, inp_model, idxl,
          //     vec_w[idx]);
          // std::size_t ridx = sem.LtG_T(idxlow_e, 0);
          // std::size_t len_ms = lentor - ridx;
          std::size_t ridx = ridxtor[idx - myff.i1()];
          std::size_t len_ms = lentor - ridx;

          // setup
          Complex w = vec_w[idx] + ieps;
          SMATRIX mat_w_tor_red =
              mat_ke_tor.block(ridx, ridx, len_ms, len_ms) -
              w * w * mat_in_tor.block(ridx, ridx, len_ms, len_ms);

          // if including attenuation
          if (params.attenuation()) {
            // std::cout << "Including attenuation in toroidal modes\n";
            mat_w_tor_red += (myi + twodivpi * std::log(vec_w[idx] / w0)) *
                             mat_ke_tor_atten.block(ridx, ridx, len_ms, len_ms);
          }
          mat_w_tor_red.makeCompressed();

          // std::cout << "Before reduced force vector\n";
          auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
          // solver1.compute(mat_w_tor_red);

          // std::cout << "Before compute\n";
          auto idxn = myff.i2() - idx - 1;
          // std::cout << "idxn: " << idxn << "\n";
          if ((idxn % nskip) == 0) {
            solver1.compute(mat_w_tor_red);
          } else {
            solver1.factorize(mat_w_tor_red);
          }
          // solver1.compute(mat_w_tor_red);

          // std::cout << "Before solve\n";
          MATRIX vec_sol = solver1.solve(f_red);
          auto lidx = lowidx - ridx;

          // std::cout << "Before compute\n";
          // compute responses
          auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 2);
          Complex mfact;
          if (params.output_type() == 0) {
            mfact = -myi / w;
          } else if (params.output_type() == 1) {
            mfact = 1.0;
          } else if (params.output_type() == 2) {
            mfact = myi * w;
          }
          vec_raw_l.col(idx) +=
              mfact * RED_C.cwiseProduct(testmult_red).rowwise().sum();
        };

#pragma omp critical(torvecadd)
        {
          vec_raw += vec_raw_l;
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // do spheroidals:
  if (inc_sph) {
    std::cout << "Doing Spheroidal Modes\n";
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;
    using namespace std::chrono;
    using clock = high_resolution_clock;

// loop over l values
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

        auto tot_setup = microseconds::zero();
        auto tot_setup_1 = microseconds::zero();
        auto tot_setup_2 = microseconds::zero();
        auto tot_setup_3 = microseconds::zero();
        auto tot_setup_4 = microseconds::zero();
        auto tot_setup_5 = microseconds::zero();
        auto tot_setup_6 = microseconds::zero();
        auto tot_setup_7 = microseconds::zero();
        // auto total_setup_1_duration = microseconds::zero();
        // auto total_setup_2_duration = microseconds::zero();
        // auto total_setup_3_duration = microseconds::zero();
        auto tot_mat = microseconds::zero();
        auto tot_factorize = microseconds::zero();
        auto tot_solve = microseconds::zero();
        auto tot_matmult = microseconds::zero();
        auto tot_fullsolve = microseconds::zero();

        timer1.start();
#pragma omp critical(sphoutput)
        {
          // std::cout << "Spheroidal,  l = " << idxl << "\n";10
        }
        auto t_ss = clock::now();
        ///////////////////////////////////////////////////////////////////////////////
        // matrices
        SMATRIX ke_s = sem.MAT_KE_S(idxl).cast<Complex>();
        SMATRIX in_s = sem.MAT_IN_S(idxl).cast<Complex>();
        SMATRIX ke_s_a = sem.MAT_KE_S_ATTEN(idxl).cast<Complex>();
        ke_s.makeCompressed();
        in_s.makeCompressed();
        ke_s_a.makeCompressed();
        auto ts1 = clock::now();
        tot_setup_1 += duration_cast<microseconds>(ts1 - t_ss);

        ///////////////////////////////////////////////////////////////////////////////
        // std::cout << "Check 1\n";

        // MATRIX RV_VALS = sem.RV_FULL(params, idxl);
        auto t_ss2 = clock::now();
        MATRIX RV_VALS = param_info.RV_FULL_SPH(idxl);
        auto t_sst2 = clock::now();
        tot_setup_2 += duration_cast<microseconds>(t_sst2 - t_ss2);

        auto t_ss3 = clock::now();
        MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);
        auto t_sst3 = clock::now();
        tot_setup_3 += duration_cast<microseconds>(t_sst3 - t_ss3);
        // auto stop_setup_2 = clock::now();
        // total_setup_2_duration +=
        //     duration_cast<microseconds>(stop_setup_2 - t_ss);
        // std::cout << "Check 3\n";
        // t_ss = clock::now();
        auto t_ss4 = clock::now();
        MATRIX RED_C = RV_VALS * F_VALS;
        auto t_sst4 = clock::now();
        tot_setup_4 += duration_cast<microseconds>(t_sst4 - t_ss4);
        // auto stop_setup_4 = clock::now();
        // auto stop_setup_3 = clock::now();
        // total_setup_3_duration +=
        //     duration_cast<microseconds>(stop_setup_3 - t_ss);
        // std::cout << "Check 4\n";
        auto t_ss5 = clock::now();
        MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
        auto t_sst5 = clock::now();
        tot_setup_5 += duration_cast<microseconds>(t_sst5 - t_ss5);

        // std::cout << "Check 5\n";
        // t_ss = clock::now();
        auto t_ss6 = clock::now();
        MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);
        auto t_sst6 = clock::now();
        tot_setup_6 += duration_cast<microseconds>(t_sst6 - t_ss6);

        // std::cout << "Check 6\n";
        auto t_ss7 = clock::now();
        // int nskip_local =
        //     nskip * static_cast<int>(std::ceil(std::pow(double(idxl),
        //     0.25)));
        int nskip_local = nskip;
        auto vec_ridx = SpectralTools::AllIndices_SPH(
            sem, inp_model, idxl, myff, idx_source, nskip_local);
        auto t_sst7 = clock::now();
        tot_setup_7 += duration_cast<microseconds>(t_sst7 - t_ss7);

        // timing
        auto t_sts = clock::now();
        tot_setup += duration_cast<microseconds>(t_sts - t_ss);

        ///////////////////////////////////////////////////////////////////////////////
        auto t_pre_mat = clock::now();

        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          auto t_sm = clock::now();
          // indices
          std::size_t idx_rs = vec_ridx[idx - myff.i1()];
          // idx_rs = 0;
          std::size_t len_ms = lensph - idx_rs;
          auto lidx = lowidx - idx_rs;

          // setup
          Complex w = vec_w[idx] + ieps;

          // Reduced matrix
          SMATRIX mat_sph = ke_s.block(idx_rs, idx_rs, len_ms, len_ms) -
                            w * w * in_s.block(idx_rs, idx_rs, len_ms, len_ms);
          if (params.attenuation()) {
            mat_sph += (myi + twodivpi * std::log(vec_w[idx] / w0)) *
                       ke_s_a.block(idx_rs, idx_rs, len_ms, len_ms);
          }
          mat_sph.makeCompressed();
          auto f_red = F_BASE.block(idx_rs, 0, len_ms, F_BASE.cols());
          auto t_stm = clock::now();
          tot_mat += duration_cast<microseconds>(t_stm - t_sm);

          /////////////////////////////////////////
          // factorize
          auto t_sf = clock::now();
          auto idxn = myff.i2() - idx - 1;
          if ((idxn % nskip_local) == 0) {
            solver1.compute(mat_sph);
          } else {
            solver1.factorize(mat_sph);
          }
          // solver1.compute(mat_sph);
          auto t_stf = clock::now();
          tot_factorize += duration_cast<microseconds>(t_stf - t_sf);

          /////////////////////////////////////////
          // solve
          MATRIX vec_sol = solver1.solve(f_red);
          auto t_stsol = clock::now();
          tot_solve += duration_cast<microseconds>(t_stsol - t_stf);

          /////////////////////////////////////////
          // compute responses
          auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);
          Complex mfact;
          if (params.output_type() == 0) {
            mfact = -myi / w;
          } else if (params.output_type() == 1) {
            mfact = 1.0;
          } else if (params.output_type() == 2) {
            mfact = myi * w;
          }
          vec_raw_l.col(idx) +=
              mfact * RED_C.cwiseProduct(testmult_red).rowwise().sum();
          auto t_stmm = clock::now();
          tot_matmult += duration_cast<microseconds>(t_stmm - t_stsol);
        };

#pragma omp critical(sphvecadd)
        {
          vec_raw += vec_raw_l;
        }
        auto t_post_mat = clock::now();
        /*
        tot_fullsolve += duration_cast<microseconds>(t_post_mat - t_pre_mat);
        std::cout << "Total time for setup: " << tot_setup.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 1: " << tot_setup_1.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 2: " << tot_setup_2.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 3: " << tot_setup_3.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 4: " << tot_setup_4.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 5: " << tot_setup_5.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 6: " << tot_setup_6.count() / 1e6
                  << " s\n";
        std::cout << "  - time for setup 7: " << tot_setup_7.count() / 1e6
                  << " s\n";
        std::cout << "Total time for solve: " << tot_fullsolve.count() / 1e6
                  << " s\n";
        // std::cout << "Total time for factorization: "
        //           << tot_factorize.count() / 1e6 << " s\n";
        // std::cout << "Total time for solving: " << tot_solve.count() / 1e6
        //           << " s\n";
        // std::cout << "Total time for matrix multiplications: "
        //           << tot_matmult.count() / 1e6 << " s\n";

        timer1.stop("Time for l = " + std::to_string(idxl));
        std::cout << "\n";
        */
      }
    }
  }

  return vec_raw;
};

template <class model1d>
auto
Sparse_F_Spec::FrequencySpectrum_Variable(SpectraSolver::FreqFull &myff,
                                          model1d &inp_model,
                                          SourceInfo::EarthquakeCMT &cmt,
                                          InputParameters &params, int NQ,
                                          int nskip, int num_chunks) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  // using SUMF = Eigen::UmfPackLU<SMATRIX>;
  using BCS = Eigen::BiCGSTAB<SMATRIX, Eigen::IncompleteLUT<Complex>>;
  using BCSP = Eigen::BiCGSTAB<SMATRIX, Eigen::FEMPreconditioner<Complex>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  // auto NQ = 6;
  // Take two elements per wavelength at the maximum frequency:
  double maxstep = 0.25 / myff.f22();
  // Full1D::specsem sem(inp_model, maxstep, NQ, params.lmax());
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  auto tref = inp_model.TREF();
  auto w0 = 2.0 * 3.14159265358979323846 / tref;
  auto twodivpi = 2.0 / 3.14159265358979323846;

  // std::cout << "Reference frequency: " << w0 << " Hz\n";
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  //////////////////////////////////////////////////////////////////////////////
  // divide the frequencies into chunks to have different maximum step sizes for
  // different frequencies. We start with
  auto wl = vec_w[myff.i1()];
  auto wh = vec_w[myff.i2() - 1];
  // std::cout << "Minimum frequency: " << wl << " Hz\n";
  // std::cout << "Maximum frequency: " << wh << " Hz\n";
  // std::cout << "Frequencies: " << myff.f11() << " to " << myff.f22() << "
  // Hz\n";
  // int num_chunks = 10;
  double powval = 1.0 / static_cast<double>(num_chunks);
  auto wrat = std::pow(wh / wl, powval);
  auto wdiff = (wh - wl) / static_cast<double>(num_chunks);
  std::vector<std::vector<double>> freq_chunks;
  std::vector<std::vector<int>> idx_chunks;
  auto idxw = myff.i1();
  for (int idx = 0; idx < num_chunks; ++idx) {
    std::vector<double> chunk;
    std::vector<int> tmp;
    // std::cout << "Frequency chunk " << idx << ": " << wl * std::pow(wrat,
    // idx)
    //           << " to " << wl * std::pow(wrat, idx + 1) << " Hz\n";
    while (vec_w[idxw] <= wl + wdiff * (idx + 1)) {
      chunk.push_back(vec_w[idxw]);
      tmp.push_back(idxw);
      ++idxw;
    }
    freq_chunks.push_back(chunk);
    idx_chunks.push_back(tmp);
  }
  for (int idx = 0; idx < num_chunks; ++idx) {
    std::cout << "Frequency chunk " << idx << ": " << freq_chunks[idx].front()
              << " to " << freq_chunks[idx].back() << " Hz\n";
  }

  // vector of maximum step sizes for each frequency chunk
  std::vector<double> max_steps;
  for (int idx = 0; idx < num_chunks; ++idx) {
    double tmp = 1.57 / freq_chunks[idx].back();
    max_steps.push_back(tmp);
  }

  // vector of sems
  std::vector<Full1D::specsem> sems;
  for (int idx = 0; idx < num_chunks; ++idx) {
    sems.emplace_back(inp_model, max_steps[idx], NQ, params.lmax());
  }

  //////////////////////////////////////////////////////////////////////////////
  timer1.start();
  SLU solver, solver1;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamInfo param_info(params, lmax);
  // auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();

  // std::cout << "Output type: " << params.output_type() << "\n";
  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    timer1.start();

    std::cout << "Doing Radial Modes\n";
    // auto idx_chunk = 0;
    for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
      Full1D::specsem &sem = sems[idx_chunk];
      auto rec_elems = sem.Receiver_Elements(params);
      auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
      auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
      int lenidx = upidx - lowidx + 1;
      SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
      SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
      SMATRIX ke_r_atten = sem.MAT_KE_R_ATTEN().cast<Complex>();
      ke_r.makeCompressed();
      in_r.makeCompressed();
      ke_r_atten.makeCompressed();

      auto twopi = 2.0 * 3.14159265358979323846;

      // calculate force vector
      MATRIX f_r = sem.CalculateForce_R(cmt);

      std::vector<MATRIX> vec_RV_Z;
      for (int idxr = 0; idxr < num_rec; ++idxr) {
        vec_RV_Z.push_back(
            sem.RV_Z_R(params, idxr).block(lowidx, 0, lenidx, 1));
      }

#pragma omp parallel default(shared) private(solver)
      {
#pragma omp for schedule(dynamic, 10)
        for (int idxw = idx_chunks[idx_chunk][0];
             idxw < idx_chunks[idx_chunk].size(); ++idxw) {
          auto wval = freq_chunks[idx_chunk][idxw];
          // complex frequency
          Complex w = wval + ieps;

          // build matrix and solve
          SMATRIX w_r = -w * w * in_r + ke_r;
          // if including attenuation
          if (params.attenuation()) {
            // std::cout << "Including attenuation in toroidal modes\n";
            w_r += (myi + twodivpi * std::log(wval / w0)) * ke_r_atten;
          }
          w_r.makeCompressed();
          solver.compute(w_r);
          MATRIX vec_x = solver.solve(f_r);

          // compute responses
          for (int idxr = 0; idxr < num_rec; ++idxr) {
            // index
            auto idxpl = 3 * idxr;
            Complex mfact;
            if (params.output_type() == 0) {
              mfact = -myi / w;
            } else if (params.output_type() == 1) {
              mfact = 1.0;
            } else if (params.output_type() == 2) {
              mfact = myi * w;
            }
// find response
#pragma omp critical(torvecadd)
            {
              vec_raw(idxpl, idxw) +=
                  mfact * vec_RV_Z[idxr]
                              .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                              .sum();
            }
          }
        };
      }
    }
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";

#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        // std::cout << "Toroidal,  l = " << idxl << "\n";
        MATRIX RV_VALS = param_info.RV_FULL_TOR(idxl);
        for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          Full1D::specsem &sem = sems[idx_chunk];
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_T(rec_elems[0], 0);
          auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
          int lenidx = upidx - lowidx;
          auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;

          // get matrices
          SMATRIX mat_ke_tor = sem.MAT_KE_T_K(idxl).cast<Complex>();
          SMATRIX mat_in_tor = sem.MAT_IN_T_K(idxl).cast<Complex>();
          SMATRIX mat_ke_tor_atten = sem.MAT_KE_T_ATTEN(idxl).cast<Complex>();
          mat_ke_tor.makeCompressed();
          mat_in_tor.makeCompressed();
          mat_ke_tor_atten.makeCompressed();
          // std::cout << "Test 1\n";
          ///////////////////////////////////////////////////////////////////////
          // setting up forces and receiver vectors
          MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);
          MATRIX RED_C = RV_VALS * F_VALS;
          MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
          MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);

          // std::cout << "Test 2\n";
          // indices of reduced matrices at each frequency
          auto ridxtor = SpectralTools::AllIndices_TOR(
              sem, inp_model, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i2 = idx_chunks[idx_chunk].back();
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l = MATRIX::Zero(3 * params.num_receivers(),
                                          freq_chunks[idx_chunk].size());

          // std::cout << "Test 3\n";

          // std::cout << "Test 4\n";

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            auto idxw = idx_chunks[idx_chunk][idx];

            // std::cout << "Frequency: " << wval << " Hz\n";

            // std::cout << "Test f.1\n";
            // lower element at this frequency
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;
            // std::cout << "Test f.2\n";
            // setup
            Complex w = wval + ieps;
            SMATRIX mat_w_tor_red =
                mat_ke_tor.block(ridx, ridx, len_ms, len_ms) -
                w * w * mat_in_tor.block(ridx, ridx, len_ms, len_ms);
            // std::cout << "Test f.3\n";
            // if including attenuation
            if (params.attenuation()) {
              mat_w_tor_red +=
                  (myi + twodivpi * std::log(wval / w0)) *
                  mat_ke_tor_atten.block(ridx, ridx, len_ms, len_ms);
            }
            mat_w_tor_red.makeCompressed();
            // std::cout << "Test f.4\n";
            auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;
            if ((idxn % nskip) == 0) {
              solver1.compute(mat_w_tor_red);
            } else {
              solver1.factorize(mat_w_tor_red);
            }
            // std::cout << "Test f.5\n";
            MATRIX vec_sol = solver1.solve(f_red);
            auto lidx = lowidx - ridx;
            // std::cout << "Test f.6\n";
            // compute responses
            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 2);
            Complex mfact;
            if (params.output_type() == 0) {
              mfact = -myi / w;
            } else if (params.output_type() == 1) {
              mfact = 1.0;
            } else if (params.output_type() == 2) {
              mfact = myi * w;
            }
            // std::cout << "Test f.7\n";
            vec_raw_l.col(idx) +=
                mfact * RED_C.cwiseProduct(testmult_red).rowwise().sum();
          };

#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), len_chunk) +=
                vec_raw_l;
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // SUMF solver2;
  BCSP solver3;
  // spheroidals
  if (inc_sph) {
    std::cout << "Doing Spheroidal Modes\n";
    using namespace std::chrono;
    using clock = high_resolution_clock;

    auto tot_rv_vals = microseconds::zero();
    auto tot_sem = microseconds::zero();
    auto tot_mat = microseconds::zero();
    auto tot_f_vals = microseconds::zero();
    auto tot_red_c = microseconds::zero();
    auto tot_f_base = microseconds::zero();
    auto tot_rv_base = microseconds::zero();
    auto tot_ridx = microseconds::zero();
    auto tot_rmat = microseconds::zero();
    auto tot_factorize = microseconds::zero();
    auto tot_solve = microseconds::zero();
    auto tot_matmult = microseconds::zero();

#pragma omp parallel default(shared) private(solver1, solver3)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        // std::cout << "Toroidal,  l = " << idxl << "\n";
        auto tot_l_rv_vals = microseconds::zero();
        auto tot_l_sem = microseconds::zero();
        auto tot_l_f_vals = microseconds::zero();
        auto tot_l_red_c = microseconds::zero();
        auto tot_l_f_base = microseconds::zero();
        auto tot_l_rv_base = microseconds::zero();
        auto tot_l_ridx = microseconds::zero();
        auto tot_l_rmat = microseconds::zero();
        auto tot_l_factorize = microseconds::zero();
        auto tot_l_solve = microseconds::zero();
        auto tot_l_matmult = microseconds::zero();

        auto t_rv_vals = clock::now();
        MATRIX RV_VALS = param_info.RV_FULL_SPH(idxl);
        auto t_rv_vals_end = clock::now();
        tot_l_rv_vals += duration_cast<microseconds>(t_rv_vals_end - t_rv_vals);

#pragma omp critical(sphoutput)
        {
          // std::cout << "Spheroidal,  l = " << idxl << "\n";
        }
        for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          // std::cout << "Frequency chunk " << idx_chunk << ": "
          //           << freq_chunks[idx_chunk].front() << " to "
          //           << freq_chunks[idx_chunk].back() << " Hz\n";
          auto t_sem = clock::now();
          Full1D::specsem &sem = sems[idx_chunk];
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
          auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          auto t_sem_end = clock::now();
          tot_l_sem += duration_cast<microseconds>(t_sem_end - t_sem);

          // get matrices
          auto t_mat = clock::now();
          SMATRIX mat_ke_sph = sem.MAT_KE_S(idxl).cast<Complex>();
          SMATRIX mat_in_sph = sem.MAT_IN_S(idxl).cast<Complex>();
          SMATRIX mat_ke_sph_atten = sem.MAT_KE_S_ATTEN(idxl).cast<Complex>();
          mat_ke_sph.makeCompressed();
          mat_in_sph.makeCompressed();
          mat_ke_sph_atten.makeCompressed();
          auto t_mat_end = clock::now();
          tot_l_rmat += duration_cast<microseconds>(t_mat_end - t_mat);
          // std::cout << "Test 1\n";
          ///////////////////////////////////////////////////////////////////////
          // setting up forces and receiver vectors
          auto t_f_vals = clock::now();
          MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);
          auto t_f_vals_end = clock::now();
          tot_l_f_vals += duration_cast<microseconds>(t_f_vals_end - t_f_vals);
          auto t_red_c = clock::now();
          MATRIX RED_C = RV_VALS * F_VALS;
          auto t_red_c_end = clock::now();
          tot_l_red_c += duration_cast<microseconds>(t_red_c_end - t_red_c);
          auto t_f_base = clock::now();
          MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
          auto t_f_base_end = clock::now();
          tot_l_f_base += duration_cast<microseconds>(t_f_base_end - t_f_base);
          auto t_rv_base = clock::now();
          MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);
          auto t_rv_base_end = clock::now();
          tot_l_rv_base +=
              duration_cast<microseconds>(t_rv_base_end - t_rv_base);
          // std::cout << "Test 2\n";
          // indices of reduced matrices at each frequency
          auto t_ridx = clock::now();
          auto ridxsph = SpectralTools::AllIndices_SPH(
              sem, inp_model, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i2 = idx_chunks[idx_chunk].back();
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l = MATRIX::Zero(3 * params.num_receivers(),
                                          freq_chunks[idx_chunk].size());
          auto t_ridx_end = clock::now();
          tot_l_ridx += duration_cast<microseconds>(t_ridx_end - t_ridx);
          // std::cout << "Test 3\n";

          // std::cout << "Test 4\n";

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            auto idxw = idx_chunks[idx_chunk][idx];

            // std::cout << "Frequency: " << wval << " Hz\n";

            // std::cout << "Test f.1\n";
            // lower element at this frequency
            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;
            // std::cout << "Test f.2\n";
            // setup
            auto t_rmat = clock::now();
            Complex w = wval + ieps;
            SMATRIX mat_w_sph_red =
                mat_ke_sph.block(ridx, ridx, len_ms, len_ms) -
                w * w * mat_in_sph.block(ridx, ridx, len_ms, len_ms);
            // std::cout << "Test f.3\n";
            // if including attenuation
            if (params.attenuation()) {
              mat_w_sph_red +=
                  (myi + twodivpi * std::log(wval / w0)) *
                  mat_ke_sph_atten.block(ridx, ridx, len_ms, len_ms);
            }
            mat_w_sph_red.makeCompressed();
            auto t_rmat_end = clock::now();
            tot_l_rmat += duration_cast<microseconds>(t_rmat_end - t_rmat);
            // std::cout << "Test f.4\n";
            auto t_f_red = clock::now();
            MATRIX f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto t_f_red_end = clock::now();
            tot_l_f_vals += duration_cast<microseconds>(t_f_red_end - t_f_red);

            auto t_factorize = clock::now();
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;

            // std::cout << "Frequency: " << wval << " Hz, idxn: " << idxn <<
            // "\n";
            // if ((idxn % nskip) == 0) {
            //   solver1.compute(mat_w_sph_red);
            // } else {
            //   solver1.factorize(mat_w_sph_red);
            // }
            // solver1.compute(mat_w_sph_red);

            solver3.preconditioner().addmatrix(mat_w_sph_red);
            Eigen::MatrixXcd vec_guess = solver3.preconditioner().solve(f_red);

            auto t_factorize_end = clock::now();
            tot_l_factorize +=
                duration_cast<microseconds>(t_factorize_end - t_factorize);
            solver3.compute(mat_w_sph_red);
            solver3.setTolerance(1e-10);
            //  std::cout << "Original
            // norm: " << vec_res.norm() << "\n";

            // std::cout << "Test f.5\n";
            auto t_solve = clock::now();
            // MATRIX vec_sol = solver1.solve(f_red);
            // MATRIX vec_sol = solver3.solve(f_red);
            MATRIX vec_sol = solver3.solveWithGuess(f_red, vec_guess);
            auto t_solve_end = clock::now();
            tot_l_solve += duration_cast<microseconds>(t_solve_end - t_solve);
            auto lidx = lowidx - ridx;
            // std::cout << "Test f.6\n";
            // compute responses
            auto t_matmult = clock::now();
            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);

            Complex mfact;
            if (params.output_type() == 0) {
              mfact = -myi / w;
            } else if (params.output_type() == 1) {
              mfact = 1.0;
            } else if (params.output_type() == 2) {
              mfact = myi * w;
            }
            // std::cout << "Test f.7\n";
            vec_raw_l.col(idx) +=
                mfact * RED_C.cwiseProduct(testmult_red).rowwise().sum();
            auto t_matmult_end = clock::now();
            tot_l_matmult +=
                duration_cast<microseconds>(t_matmult_end - t_matmult);
          };

#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), len_chunk) +=
                vec_raw_l;
          }
        }
#pragma omp critical(sphoutput)
        {
          tot_rv_vals += tot_l_rv_vals;
          tot_sem += tot_l_sem;
          tot_f_vals += tot_l_f_vals;
          tot_red_c += tot_l_red_c;
          tot_f_base += tot_l_f_base;
          tot_rv_base += tot_l_rv_base;
          tot_ridx += tot_l_ridx;
          tot_rmat += tot_l_rmat;
          tot_factorize += tot_l_factorize;
          tot_solve += tot_l_solve;
          tot_matmult += tot_l_matmult;
          // std::cout << "Spheroidal,  l = " << idxl << "\n";
        }
      }
    }

    // output timing information
    std::cout << "Total time for RV_VALS: " << tot_rv_vals.count() / 1e6
              << " s\n";
    std::cout << "Total time for SEM construction: " << tot_sem.count() / 1e6
              << " s\n";
    std::cout << "Total time for F_VALS: " << tot_f_vals.count() / 1e6
              << " s\n";
    std::cout << "Total time for RED_C: " << tot_red_c.count() / 1e6 << " s\n";
    std::cout << "Total time for F_BASE: " << tot_f_base.count() / 1e6
              << " s\n";
    std::cout << "Total time for RV_BASE: " << tot_rv_base.count() / 1e6
              << " s\n";
    std::cout << "Total time for RIDX: " << tot_ridx.count() / 1e6 << " s\n";
    std::cout << "Total time for reduced matrix construction: "
              << tot_rmat.count() / 1e6 << " s\n";
    std::cout << "Total time for factorization: " << tot_factorize.count() / 1e6
              << " s\n";
    std::cout << "Total time for solve: " << tot_solve.count() / 1e6 << " s\n";
    std::cout << "Total time for matrix multiplication: "
              << tot_matmult.count() / 1e6 << " s\n";
  }

  return vec_raw;
};   // namespace SPARSESPEC

template <class model1d>
auto
Sparse_F_Spec::FrequencySpectrum_RED(SpectraSolver::FreqFull &myff,
                                     model1d &inp_model,
                                     SourceInfo::EarthquakeCMT &cmt,
                                     InputParameters &params, int NQ, int nskip,
                                     int num_chunks, SRInfo &srinfo) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  // using SUMF = Eigen::UmfPackLU<SMATRIX>;
  using BCS = Eigen::BiCGSTAB<SMATRIX, Eigen::IncompleteLUT<Complex>>;
  using BCSP = Eigen::BiCGSTAB<SMATRIX, Eigen::FEMPreconditioner<Complex>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  // auto NQ = 6;
  // Take two elements per wavelength at the maximum frequency:
  double maxstep = 0.25 / myff.f22();
  // Full1D::specsem sem(inp_model, maxstep, NQ, params.lmax());
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  auto tref = inp_model.TREF();
  auto w0 = 2.0 * 3.14159265358979323846 / tref;
  auto twodivpi = 2.0 / 3.14159265358979323846;

  // std::cout << "Reference frequency: " << w0 << " Hz\n";
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  //////////////////////////////////////////////////////////////////////////////
  // divide the frequencies into chunks to have different maximum step sizes for
  // different frequencies. We start with
  auto wl = vec_w[myff.i1()];
  auto wh = vec_w[myff.i2() - 1];
  // std::cout << "Minimum frequency: " << wl << " Hz\n";
  // std::cout << "Maximum frequency: " << wh << " Hz\n";
  // std::cout << "Frequencies: " << myff.f11() << " to " << myff.f22() << "
  // Hz\n";
  // int num_chunks = 10;
  double powval = 1.0 / static_cast<double>(num_chunks);
  auto wrat = std::pow(wh / wl, powval);
  auto wdiff = (wh - wl) / static_cast<double>(num_chunks);
  std::vector<std::vector<double>> freq_chunks;
  std::vector<std::vector<int>> idx_chunks;
  auto idxw = myff.i1();
  for (int idx = 0; idx < num_chunks; ++idx) {
    std::vector<double> chunk;
    std::vector<int> tmp;
    // std::cout << "Frequency chunk " << idx << ": " << wl * std::pow(wrat,
    // idx)
    //           << " to " << wl * std::pow(wrat, idx + 1) << " Hz\n";
    while (vec_w[idxw] <= wl + wdiff * (idx + 1)) {
      chunk.push_back(vec_w[idxw]);
      tmp.push_back(idxw);
      ++idxw;
    }
    freq_chunks.push_back(chunk);
    idx_chunks.push_back(tmp);
  }
  for (int idx = 0; idx < num_chunks; ++idx) {
    std::cout << "Frequency chunk " << idx << ": " << freq_chunks[idx].front()
              << " to " << freq_chunks[idx].back() << " Hz\n";
  }

  // vector of maximum step sizes for each frequency chunk
  std::vector<double> max_steps;
  for (int idx = 0; idx < num_chunks; ++idx) {
    double tmp = 1.57 / freq_chunks[idx].back();
    max_steps.push_back(tmp);
  }

  // vector of sems
  std::vector<Full1D::specsem> sems;
  for (int idx = 0; idx < num_chunks; ++idx) {
    sems.emplace_back(inp_model, max_steps[idx], NQ, params.lmax());
  }

  //////////////////////////////////////////////////////////////////////////////
  timer1.start();
  SLU solver, solver1;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamRedInfo param_info(params, lmax);
  // auto NQ = sem.mesh().NN();

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
  // std::cout << "mtype: " << mtype << "\n";
  if (mtype == 4) {
    inc_rad = true;
    inc_tor = true;
    inc_sph = true;
  } else if (mtype == 1) {
    inc_rad = true;
  } else if (mtype == 2) {
    inc_tor = true;
  } else if (mtype == 3) {
    inc_sph = true;
  }
  if (lmin > 0) {
    inc_rad = false;
  }
  if (lmax < 1) {
    inc_tor = false;
    inc_sph = false;
  }

  // change lmin
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();

  // std::cout << "Output type: " << params.output_type() << "\n";
  ///////////////////////////////////

  //   timer1.start();

  // do radials
  if (inc_rad) {
    timer1.start();

    std::cout << "Doing Radial Modes\n";
    // auto idx_chunk = 0;
    for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
      // std::cout << "Frequency chunk " << idx_chunk << ": "
      //           << freq_chunks[idx_chunk].front() << " to "
      //           << freq_chunks[idx_chunk].back() << " Hz\n\n";
      // std::cout << "idxw from: " << idx_chunks[idx_chunk][0] << " to "
      //           << idx_chunks[idx_chunk].back() << "\n";
      Full1D::specsem &sem = sems[idx_chunk];
      auto rec_elems = sem.Receiver_Elements(params);
      auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
      auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
      int lenidx = upidx - lowidx + 1;

      // std::cout << "Lowidx: " << lowidx << ", Upidx: " << upidx
      //           << ", Lenidx: " << lenidx << "\n";
      SMATRIX ke_r = sem.MAT_KE_R().cast<Complex>();
      SMATRIX in_r = sem.MAT_IN_R().cast<Complex>();
      SMATRIX ke_r_atten = sem.MAT_KE_R_ATTEN().cast<Complex>();
      ke_r.makeCompressed();
      in_r.makeCompressed();
      ke_r_atten.makeCompressed();

      auto twopi = 2.0 * 3.14159265358979323846;

      // calculate force vector
      // MATRIX f_r = sem.CalculateForce_Red_R(cmt);

      // calculate force vector
      MATRIX f_r = sem.CalculateForce_Red_R(cmt);

      MATRIX vec_red_z = sem.RV_RED_Z_R(params);
      // std::cout << "idx_chunk: " << idx_chunk << "\n";
      // std::cout << vec_red_z << "\n\n";
      // auto testrec = sem.RV_Z_R(params, 0).block(lowidx, 0, lenidx, 1);
      // std::cout << "Test rec: \n" << testrec << "\n\n";
      // #pragma omp parallel default(shared) private(solver)
      {
        // #pragma omp for schedule(dynamic, 10)
        for (int idx = 0; idx < idx_chunks[idx_chunk].size(); ++idx) {
          double wval = freq_chunks[idx_chunk][idx];
          int idxw = idx_chunks[idx_chunk][idx];
          // std::cout << "idx_chunk: " << idx_chunk << ", frequency: " << wval
          //           << " Hz\n";
          // complex frequency
          Complex w = wval + ieps;

          // build matrix and solve
          SMATRIX w_r = -w * w * in_r + ke_r;
          // if including attenuation
          if (params.attenuation()) {
            // std::cout << "Including attenuation in toroidal modes\n";
            w_r += (myi + twodivpi * std::log(wval / w0)) * ke_r_atten;
          }

          w_r.makeCompressed();
          solver.compute(w_r);
          MATRIX vec_x = solver.solve(f_r);
          auto tmp = vec_x.block(lowidx, 0, lenidx, 1);

          Complex mfact;
          if (params.output_type() == 0) {
            mfact = -myi / w;
          } else if (params.output_type() == 1) {
            mfact = 1.0;
          } else if (params.output_type() == 2) {
            mfact = myi * w;
          }
          auto resval = mfact * (vec_red_z.transpose() * tmp).sum();
          // if ((wval > 4.0) && (wval < 7.0)) {
          //   std::cout << "Test solution at frequency: " << wval << " Hz\n";
          //   std::cout << std::setprecision(4) << tmp.transpose() << "\n";
          //   std::cout << std::setprecision(4)
          //             << "Test reduced receiver vector: \n"
          //             << vec_red_z.transpose() << "\n";
          //   std::cout << std::setprecision(4)
          //             << "Test response value: " << resval << "\n";
          // }
          // compute responses
          for (int idxr = 0; idxr < num_rec; ++idxr) {
            // index
            auto idxpl = 3 * idxr;
            // if ((wval > 4.0) && (wval < 7.0)) {
            //   std::cout << "Adding response to vec_raw at index: (" << idxpl
            //             << ", " << idxw << ")\n";
            // }
            // find response
            // #pragma omp critical(torvecadd)
            {
              vec_raw(idxpl, idxw) += resval;
            }
          }
          // if ((wval > 4.0) && (wval < 7.0)) {
          //   std::cout << "\n";
          // }
        };
      }
    }
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";

#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        // std::cout << "Toroidal,  l = " << idxl << "\n";
        MATRIX RV_VALS = srinfo.RV_RED_TOR(idxl);
        for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          Full1D::specsem &sem = sems[idx_chunk];
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_T(rec_elems[0], 0);
          auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
          int lenidx = upidx - lowidx;
          auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;

          // get matrices
          SMATRIX mat_ke_tor = sem.MAT_KE_T_K(idxl).cast<Complex>();
          SMATRIX mat_in_tor = sem.MAT_IN_T_K(idxl).cast<Complex>();
          SMATRIX mat_ke_tor_atten = sem.MAT_KE_T_ATTEN(idxl).cast<Complex>();
          mat_ke_tor.makeCompressed();
          mat_in_tor.makeCompressed();
          mat_ke_tor_atten.makeCompressed();
          // std::cout << "Test 1\n";
          ///////////////////////////////////////////////////////////////////////
          // setting up forces and receiver vectors
          MATRIX F_VALS = sem.CalculateForce_RED_Coefficients_T(cmt, idxl, 0.0);
          MATRIX RED_C = RV_VALS * F_VALS;
          MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
          MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);

          // std::cout << "Test 2\n";
          // indices of reduced matrices at each frequency
          auto ridxtor = SpectralTools::AllIndices_TOR(
              sem, inp_model, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i2 = idx_chunks[idx_chunk].back();
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l = MATRIX::Zero(3 * params.num_receivers(),
                                          freq_chunks[idx_chunk].size());

          // std::cout << "Test 3\n";

          // std::cout << "Test 4\n";

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            auto idxw = idx_chunks[idx_chunk][idx];

            // std::cout << "Frequency: " << wval << " Hz\n";

            // std::cout << "Test f.1\n";
            // lower element at this frequency
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;
            // std::cout << "Test f.2\n";
            // setup
            Complex w = wval + ieps;
            SMATRIX mat_w_tor_red =
                mat_ke_tor.block(ridx, ridx, len_ms, len_ms) -
                w * w * mat_in_tor.block(ridx, ridx, len_ms, len_ms);
            // std::cout << "Test f.3\n";
            // if including attenuation
            if (params.attenuation()) {
              mat_w_tor_red +=
                  (myi + twodivpi * std::log(wval / w0)) *
                  mat_ke_tor_atten.block(ridx, ridx, len_ms, len_ms);
            }
            mat_w_tor_red.makeCompressed();
            // std::cout << "Test f.4\n";
            auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;
            if ((idxn % nskip) == 0) {
              solver1.compute(mat_w_tor_red);
            } else {
              solver1.factorize(mat_w_tor_red);
            }
            // std::cout << "Test f.5\n";
            MATRIX vec_sol = solver1.solve(f_red);
            auto lidx = lowidx - ridx;
            // std::cout << "Test f.6\n";
            // compute responses
            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 2);
            Complex mfact;
            if (params.output_type() == 0) {
              mfact = -myi / w;
            } else if (params.output_type() == 1) {
              mfact = 1.0;
            } else if (params.output_type() == 2) {
              mfact = myi * w;
            }
            // std::cout << "Test f.7\n";
            vec_raw_l.col(idx) +=
                mfact * RED_C.cwiseProduct(testmult_red).rowwise().sum();
          };

#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), len_chunk) +=
                vec_raw_l;
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // SUMF solver2;
  // BCSP solver3;
  // spheroidals
  if (inc_sph) {
    std::cout << "\nDoing Spheroidal Modes\n";
    using namespace std::chrono;
    using clock = high_resolution_clock;

    auto tot_rv_vals = microseconds::zero();
    auto tot_sem = microseconds::zero();
    auto tot_mat = microseconds::zero();
    auto tot_f_vals = microseconds::zero();
    auto tot_red_c = microseconds::zero();
    auto tot_f_base = microseconds::zero();
    auto tot_rv_base = microseconds::zero();
    auto tot_ridx = microseconds::zero();
    auto tot_rmat = microseconds::zero();
    auto tot_factorize = microseconds::zero();
    auto tot_solve = microseconds::zero();
    auto tot_matmult = microseconds::zero();

#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        // std::cout << "Toroidal,  l = " << idxl << "\n";
        auto tot_l_rv_vals = microseconds::zero();
        auto tot_l_sem = microseconds::zero();
        auto tot_l_f_vals = microseconds::zero();
        auto tot_l_red_c = microseconds::zero();
        auto tot_l_f_base = microseconds::zero();
        auto tot_l_rv_base = microseconds::zero();
        auto tot_l_ridx = microseconds::zero();
        auto tot_l_rmat = microseconds::zero();
        auto tot_l_factorize = microseconds::zero();
        auto tot_l_solve = microseconds::zero();
        auto tot_l_matmult = microseconds::zero();

        // #pragma omp critical(sphoutput)
        //         {
        //           std::cout << "Spheroidal,  l = " << idxl << "\n";
        //           std::cout << "Before RV_VALS\n\n";
        //         }

        // std::cout << "Test 0\n";
        auto t_rv_vals = clock::now();
        MATRIX RV_VALS = srinfo.RV_RED_SPH(idxl);   // receiver
        auto t_rv_vals_end = clock::now();
        tot_l_rv_vals += duration_cast<microseconds>(t_rv_vals_end - t_rv_vals);
        // std::cout << "Test 0.5\n";
        for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          // std::cout << "Frequency chunk " << idx_chunk << ": "
          //           << freq_chunks[idx_chunk].front() << " to "
          //           << freq_chunks[idx_chunk].back() << " Hz\n";
          auto t_sem = clock::now();
          Full1D::specsem &sem = sems[idx_chunk];
          // std::cout << "Test 0.6\n";
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
          auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          auto t_sem_end = clock::now();
          tot_l_sem += duration_cast<microseconds>(t_sem_end - t_sem);
          // #pragma omp critical(sphoutput)
          //           {
          //             std::cout << "Test 0\n";
          //           }
          // get matrices
          auto t_mat = clock::now();
          SMATRIX mat_ke_sph = sem.MAT_KE_S(idxl).cast<Complex>();
          SMATRIX mat_in_sph = sem.MAT_IN_S(idxl).cast<Complex>();
          SMATRIX mat_ke_sph_atten = sem.MAT_KE_S_ATTEN(idxl).cast<Complex>();
          mat_ke_sph.makeCompressed();
          mat_in_sph.makeCompressed();
          mat_ke_sph_atten.makeCompressed();
          auto t_mat_end = clock::now();
          tot_l_rmat += duration_cast<microseconds>(t_mat_end - t_mat);
          // std::cout << "Test 1\n";
          ///////////////////////////////////////////////////////////////////////

          // #pragma omp critical(sphoutput)
          //           {
          //             std::cout << "Test 1\n";
          //           }
          // setting up forces and receiver vectors
          auto t_f_vals = clock::now();
          MATRIX F_VALS =
              sem.CalculateForce_RED_Coefficients(cmt, idxl, 0.0);   // source
          // #pragma omp critical(sphoutput)
          //           {
          //             std::cout << "Spheroidal,  l = " << idxl << "\n";
          //             std::cout << "Size of F_VALS: " << F_VALS.rows() << " x
          //             "
          //                       << F_VALS.cols() << "\n";
          //             std::cout << "Size of RV_VALS: " << RV_VALS.rows() << "
          //             x "
          //                       << RV_VALS.cols() << "\n";
          //           }
          auto t_f_vals_end = clock::now();
          tot_l_f_vals += duration_cast<microseconds>(t_f_vals_end - t_f_vals);
          auto t_red_c = clock::now();
          MATRIX RED_C = RV_VALS * F_VALS;
          auto t_red_c_end = clock::now();
          tot_l_red_c += duration_cast<microseconds>(t_red_c_end - t_red_c);
          auto t_f_base = clock::now();
          MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
          auto t_f_base_end = clock::now();
          tot_l_f_base += duration_cast<microseconds>(t_f_base_end - t_f_base);
          auto t_rv_base = clock::now();
          MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);
          auto t_rv_base_end = clock::now();
          tot_l_rv_base +=
              duration_cast<microseconds>(t_rv_base_end - t_rv_base);
          // std::cout << "Test 2\n";
          // indices of reduced matrices at each frequency
          auto t_ridx = clock::now();
          auto ridxsph = SpectralTools::AllIndices_SPH(
              sem, inp_model, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i2 = idx_chunks[idx_chunk].back();
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l = MATRIX::Zero(3 * params.num_receivers(),
                                          freq_chunks[idx_chunk].size());
          auto t_ridx_end = clock::now();
          tot_l_ridx += duration_cast<microseconds>(t_ridx_end - t_ridx);
          // std::cout << "Test 3\n";

          // std::cout << "Test 4\n";

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            auto idxw = idx_chunks[idx_chunk][idx];

            // std::cout << "Frequency: " << wval << " Hz\n";

            // std::cout << "Test f.1\n";
            // lower element at this frequency
            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;
            // std::cout << "Test f.2\n";
            // setup
            auto t_rmat = clock::now();
            Complex w = wval + ieps;
            SMATRIX mat_w_sph_red =
                mat_ke_sph.block(ridx, ridx, len_ms, len_ms) -
                w * w * mat_in_sph.block(ridx, ridx, len_ms, len_ms);
            // std::cout << "Test f.3\n";
            // if including attenuation
            if (params.attenuation()) {
              mat_w_sph_red +=
                  (myi + twodivpi * std::log(wval / w0)) *
                  mat_ke_sph_atten.block(ridx, ridx, len_ms, len_ms);
            }
            mat_w_sph_red.makeCompressed();
            auto t_rmat_end = clock::now();
            tot_l_rmat += duration_cast<microseconds>(t_rmat_end - t_rmat);
            // std::cout << "Test f.4\n";
            auto t_f_red = clock::now();
            MATRIX f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto t_f_red_end = clock::now();
            tot_l_f_vals += duration_cast<microseconds>(t_f_red_end - t_f_red);

            auto t_factorize = clock::now();
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;

            // std::cout << "Frequency: " << wval << " Hz, idxn: " << idxn <<
            // "\n";
            if ((idxn % nskip) == 0) {
              solver1.compute(mat_w_sph_red);
            } else {
              solver1.factorize(mat_w_sph_red);
            }
            // solver1.compute(mat_w_sph_red);

            // solver3.preconditioner().addmatrix(mat_w_sph_red);
            // Eigen::MatrixXcd vec_guess =
            // solver3.preconditioner().solve(f_red);

            auto t_factorize_end = clock::now();
            tot_l_factorize +=
                duration_cast<microseconds>(t_factorize_end - t_factorize);
            // solver3.compute(mat_w_sph_red);
            // solver3.setTolerance(1e-10);
            //  std::cout << "Original
            // norm: " << vec_res.norm() << "\n";

            // std::cout << "Test f.5\n";
            auto t_solve = clock::now();
            MATRIX vec_sol = solver1.solve(f_red);
            // MATRIX vec_sol = solver3.solve(f_red);
            // MATRIX vec_sol = solver3.solveWithGuess(f_red, vec_guess);
            auto t_solve_end = clock::now();
            tot_l_solve += duration_cast<microseconds>(t_solve_end - t_solve);
            auto lidx = lowidx - ridx;
            // std::cout << "Test f.6\n";
            // compute responses
            auto t_matmult = clock::now();
            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);

            Complex mfact;
            if (params.output_type() == 0) {
              mfact = -myi / w;
            } else if (params.output_type() == 1) {
              mfact = 1.0;
            } else if (params.output_type() == 2) {
              mfact = myi * w;
            }
            // std::cout << "Test f.7\n";
            vec_raw_l.col(idx) +=
                mfact * RED_C.cwiseProduct(testmult_red).rowwise().sum();
            auto t_matmult_end = clock::now();
            tot_l_matmult +=
                duration_cast<microseconds>(t_matmult_end - t_matmult);
          };

#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), len_chunk) +=
                vec_raw_l;
          }
        }
#pragma omp critical(sphoutput)
        {
          tot_rv_vals += tot_l_rv_vals;
          tot_sem += tot_l_sem;
          tot_f_vals += tot_l_f_vals;
          tot_red_c += tot_l_red_c;
          tot_f_base += tot_l_f_base;
          tot_rv_base += tot_l_rv_base;
          tot_ridx += tot_l_ridx;
          tot_rmat += tot_l_rmat;
          tot_factorize += tot_l_factorize;
          tot_solve += tot_l_solve;
          tot_matmult += tot_l_matmult;
          // std::cout << "Spheroidal,  l = " << idxl << "\n";
        }
      }
    }

    // output timing information
    std::cout << "Total time for RV_VALS: " << tot_rv_vals.count() / 1e6
              << " s\n";
    std::cout << "Total time for SEM construction: " << tot_sem.count() / 1e6
              << " s\n";
    std::cout << "Total time for F_VALS: " << tot_f_vals.count() / 1e6
              << " s\n";
    std::cout << "Total time for RED_C: " << tot_red_c.count() / 1e6 << " s\n";
    std::cout << "Total time for F_BASE: " << tot_f_base.count() / 1e6
              << " s\n";
    std::cout << "Total time for RV_BASE: " << tot_rv_base.count() / 1e6
              << " s\n";
    std::cout << "Total time for RIDX: " << tot_ridx.count() / 1e6 << " s\n";
    std::cout << "Total time for reduced matrix construction: "
              << tot_rmat.count() / 1e6 << " s\n";
    std::cout << "Total time for factorization: " << tot_factorize.count() / 1e6
              << " s\n";
    std::cout << "Total time for solve: " << tot_solve.count() / 1e6 << " s\n";
    std::cout << "Total time for matrix multiplication: "
              << tot_matmult.count() / 1e6 << " s\n";
  }

  // do rotation:
  for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
    auto baz = srinfo.backazimuths(idxr);
    auto sbaz = std::sin(baz);
    auto cbaz = std::cos(baz);
    for (int j = 0; j < vec_raw.cols(); ++j) {
      auto tmp1 = vec_raw(3 * idxr + 1, j);
      auto tmp2 = vec_raw(3 * idxr + 2, j);
      vec_raw(3 * idxr + 1, j) = -tmp1 * cbaz + tmp2 * sbaz;
      vec_raw(3 * idxr + 2, j) = tmp1 * sbaz + tmp2 * cbaz;
    }
  }
  return vec_raw;
};   // namespace SPARSESPEC
};   // namespace SPARSESPEC

// SMATRIX mat_w_sph = ke_s;
// mat_w_sph -= w * w * in_s;
// SMATRIX mat_w_sph(ke_s.rows(), ke_s.cols());

// std::cout << vec_ke_idx.size() << ", " << vec_in_idx.size() << "\n";
// std::cout << vec_ke_idx[0][0] << ", " << vec_ke_idx[0][1] << ", "
//           << vec_ke_val[0] << "\n";
// std::cout << ke_s.coeff(vec_ke_idx[0][0], vec_ke_idx[0][1])
//           << "\n";

// for (std::size_t k = 0; k < vec_in_idx.size(); ++k) {
//   tpl_w[numke + k] =
//       T(vec_in_idx[k][0], vec_in_idx[k][1], -w * w * vec_in_val[k]);
// }
// mat_w_sph.resize(ke_s.rows(), ke_s.cols());
// mat_w_sph.setFromTriplets(tpl_w.begin(), tpl_w.end());
// mat_check.makeCompressed();
// for (int i = 0; i < 5; ++i) {
//   for (int j = 0; j < 5; ++j) {
//     std::cout << "mat_w_sph(" << i << "," << j
//               << "): " << mat_w_sph.coeff(i, j) << ", mat_check(" <<
//               i
//               << "," << j << "): " << mat_check.coeff(i, j) << "\n";
//   }
// }
// std::cout << "Difference in norms: " << (mat_w_sph -
// mat_check).norm()
//           << "\n";
// std::cout << "\nidxl: " << idxl << ", freq idx: " << idx << "\n";