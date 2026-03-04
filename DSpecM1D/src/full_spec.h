#ifndef FULL_SPEC_GUARD_H
#define FULL_SPEC_GUARD_H

#include <iostream>
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include "read_station.h"
#include "SourceInfo.h"
#include "start_element.h"
#include "spectra_master.h"
#include "ParamInfo.h"
#include <omp.h>
#include "FEM_Preconditioner.h"
#include "BiCGSTABT.h"
#include "ParamRedInfo.h"
#include "SR_Info.h"

namespace SPARSESPEC {

class Sparse_F_Spec {
public:
  Sparse_F_Spec() {};
  ~Sparse_F_Spec() {};

  template <class model1d>
  auto Spectra(SpectraSolver::FreqFull &, Full1D::specsem &, model1d &,
               SourceInfo::EarthquakeCMT &, InputParameters &, int = 10);

  template <class model1d>
  auto Spectra(SpectraSolver::FreqFull &, model1d &,
               SourceInfo::EarthquakeCMT &, InputParameters &, int, SRInfo &,
               double = 1e-4);

private:
};

template <class model1d>
auto
Sparse_F_Spec::Spectra(SpectraSolver::FreqFull &myff, Full1D::specsem &sem,
                       model1d &inp_model, SourceInfo::EarthquakeCMT &cmt,
                       InputParameters &params, int nskip) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  auto tref = inp_model.TREF();
  auto twopid = 2.0 * 3.14159265358979323846;
  auto w0 = twopid / tref;
  auto twodivpi = 2.0 / 3.14159265358979323846;

  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  SLU solver, solver1;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamInfo param_info(params, lmax);
  auto NQ = sem.mesh().NN();

  // decide which modes to include
  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
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

  // change lmin for toroidal and spheroidal modes
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();

  ///////////////////////////////////
  auto idx_source = sem.Source_Element(cmt);
  //   timer1.start();

  // do radials
  if (inc_rad) {
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
    auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SMATRIX ke_r = sem.H_R().cast<Complex>();
    SMATRIX in_r = sem.P_R().cast<Complex>();
    SMATRIX ke_r_atten = sem.H_RA().cast<Complex>();
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
          Complex mfact = 1.0;
          // if (params.output_type() == 0) {
          //   mfact = -myi / w;
          // } else if (params.output_type() == 1) {
          //   mfact = 1.0;
          // } else if (params.output_type() == 2) {
          //   mfact = myi * w;
          // }
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
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // toroidals
  if (inc_tor) {
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_T(rec_elems[0], 0);
    auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;
    auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

        // get matrices
        SMATRIX H_tor = sem.H_TK(idxl).cast<Complex>();
        SMATRIX P_tor = sem.P_TK(idxl).cast<Complex>();
        SMATRIX H_tor_atten = sem.H_TA(idxl).cast<Complex>();
        H_tor.makeCompressed();
        P_tor.makeCompressed();
        H_tor_atten.makeCompressed();

        ///////////////////////////////////////////////////////////////////////////
        // setting up forces and receiver vectors
        MATRIX RV_VALS = param_info.RV_FULL_TOR(idxl);
        MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);
        MATRIX RED_C = RV_VALS * F_VALS;
        MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
        MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);

        // indices of reduced matrices at each frequency
        auto ridxtor = SpectralTools::AllIndices_TOR(sem, inp_model, idxl, myff,
                                                     idx_source, nskip);
        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {

          // lower element at this frequency
          std::size_t ridx = ridxtor[idx - myff.i1()];
          std::size_t len_ms = lentor - ridx;

          // setup
          Complex w = vec_w[idx] + ieps;
          SMATRIX mat_w_tor_red =
              H_tor.block(ridx, ridx, len_ms, len_ms) -
              w * w * P_tor.block(ridx, ridx, len_ms, len_ms);

          // if including attenuation
          if (params.attenuation()) {
            mat_w_tor_red += (myi + twodivpi * std::log(vec_w[idx] / w0)) *
                             H_tor_atten.block(ridx, ridx, len_ms, len_ms);
          }
          mat_w_tor_red.makeCompressed();

          auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
          auto idxn = myff.i2() - idx - 1;
          // std::cout << "idxn: " << idxn << "\n";
          if ((idxn % nskip) == 0) {
            solver1.compute(mat_w_tor_red);
          } else {
            solver1.factorize(mat_w_tor_red);
          }

          MATRIX vec_sol = solver1.solve(f_red);
          auto lidx = lowidx - ridx;

          // compute responses
          auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 2);
          vec_raw_l.col(idx) +=
              RED_C.cwiseProduct(testmult_red).rowwise().sum();
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
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;

// loop over l values
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

        ///////////////////////////////////////////////////////////////////////////////
        // matrices
        SMATRIX h_s = sem.H_S(idxl).cast<Complex>();
        SMATRIX p_s = sem.P_S(idxl).cast<Complex>();
        SMATRIX h_sa = sem.H_SA(idxl).cast<Complex>();
        h_s.makeCompressed();
        p_s.makeCompressed();
        h_sa.makeCompressed();

        ///////////////////////////////////////////////////////////////////////////////
        MATRIX RV_VALS = param_info.RV_FULL_SPH(idxl);
        MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);
        MATRIX RED_C = RV_VALS * F_VALS;
        MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
        MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);

        auto vec_ridx = SpectralTools::AllIndices_SPH(sem, inp_model, idxl,
                                                      myff, idx_source, nskip);

        ///////////////////////////////////////////////////////////////////////////////
        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          // indices
          std::size_t idx_rs = vec_ridx[idx - myff.i1()];
          // idx_rs = 0;
          std::size_t len_ms = lensph - idx_rs;
          auto lidx = lowidx - idx_rs;

          // setup
          Complex w = vec_w[idx] + ieps;

          // Reduced matrix
          SMATRIX mat_sph = h_s.block(idx_rs, idx_rs, len_ms, len_ms) -
                            w * w * p_s.block(idx_rs, idx_rs, len_ms, len_ms);
          if (params.attenuation()) {
            mat_sph += (myi + twodivpi * std::log(vec_w[idx] / w0)) *
                       h_sa.block(idx_rs, idx_rs, len_ms, len_ms);
          }
          mat_sph.makeCompressed();
          auto f_red = F_BASE.block(idx_rs, 0, len_ms, F_BASE.cols());

          /////////////////////////////////////////
          // factorize
          auto idxn = myff.i2() - idx - 1;
          if ((idxn % nskip) == 0) {
            solver1.compute(mat_sph);
          } else {
            solver1.factorize(mat_sph);
          }

          /////////////////////////////////////////
          // solve
          MATRIX vec_sol = solver1.solve(f_red);

          /////////////////////////////////////////
          // compute responses
          auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);
          vec_raw_l.col(idx) +=
              RED_C.cwiseProduct(testmult_red).rowwise().sum();
        };

#pragma omp critical(sphvecadd)
        {
          vec_raw += vec_raw_l;
        }
      }
    }
  }

  // multiply by appropriate factor for frequency
  for (int idx = 0; idx < vec_raw.cols(); ++idx) {
    Complex mfact;
    if (params.output_type() == 0) {
      mfact = -myi / vec_w[idx];
    } else if (params.output_type() == 1) {
      mfact = 1.0;
    } else if (params.output_type() == 2) {
      mfact = myi * vec_w[idx];
    }
    vec_raw.col(idx) *= mfact;
  }
  return vec_raw;
};

template <class model1d>
auto
Sparse_F_Spec::Spectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                       SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                       int NQ, SRInfo &srinfo, double relerr) {
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
  auto vec_w = myff.w();
  Complex myi = Complex(0.0, 1.0);
  Complex ieps = -myff.ep() * myi;
  auto tref = inp_model.TREF();
  auto twopid = 2.0 * 3.14159265358979323846;
  auto w0 = twopid / tref;
  auto twodivpi = 2.0 / 3.14159265358979323846;
  auto num_rec = params.num_receivers();
  int nskip = (myff.i2() - myff.i1()) / 20;
  int num_chunks = std::floor((myff.f22() - myff.f11()) / 10.0) + 1;

  // response
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  //////////////////////////////////////////////////////////////////////////////
  // divide the frequencies into chunks to have different maximum step sizes for
  // different frequencies. We start with
  auto wl = vec_w[myff.i1()];
  auto wh = vec_w[myff.i2() - 1];
  auto wdiff = (wh - wl) / static_cast<double>(num_chunks);
  std::vector<std::vector<double>> freq_chunks;
  std::vector<std::vector<int>> idx_chunks;
  auto idxw = myff.i1();
  for (int idx = 0; idx < num_chunks; ++idx) {
    std::vector<double> chunk;
    std::vector<int> tmp;
    while (vec_w[idxw] <= wl + wdiff * (idx + 1)) {
      chunk.push_back(vec_w[idxw]);
      tmp.push_back(idxw);
      ++idxw;
    }
    freq_chunks.push_back(chunk);
    idx_chunks.push_back(tmp);
  }

  // vector of maximum step sizes for each frequency chunk
  // note max step size is based on the highest frequency in the chunk
  // we set maximum step size as 0.1, which is smaller than will be given for
  // the lowest frequency generally
  std::vector<double> max_steps;
  double base_len = 0.78 * 1000 / inp_model.TimeNorm() * twopid;
  double newlen = std::pow(100.0 * relerr, 1.0 / (NQ - 1.0)) * base_len;
  for (int idx = 0; idx < num_chunks; ++idx) {
    double tmp = newlen / freq_chunks[idx].back();
    if (tmp > 0.05) {
      tmp = 0.05;
    }
    max_steps.push_back(tmp);
  }
  for (auto idx : max_steps) {
    std::cout << "Max step: " << idx << "\n";
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

  bool inc_rad = false, inc_tor = false, inc_sph = false;
  auto mtype = params.type();
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

  // do radials
  if (inc_rad) {
    timer1.start();
    for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
      Full1D::specsem &sem = sems[idx_chunk];
      auto rec_elems = sem.Receiver_Elements(params);
      auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
      auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
      int lenidx = upidx - lowidx + 1;

      SMATRIX ke_r = sem.H_R().cast<Complex>();
      SMATRIX in_r = sem.P_R().cast<Complex>();
      SMATRIX ke_r_atten = sem.H_RA().cast<Complex>();
      ke_r.makeCompressed();
      in_r.makeCompressed();
      ke_r_atten.makeCompressed();

      // calculate force vector and receiver vector
      MATRIX f_r = sem.CalculateForce_Red_R(cmt);
      MATRIX vec_red_z = sem.RV_RED_Z_R(params);

#pragma omp parallel default(shared) private(solver)
      {
#pragma omp for schedule(dynamic)
        for (int idx = 0; idx < idx_chunks[idx_chunk].size(); ++idx) {
          double wval = freq_chunks[idx_chunk][idx];
          int idxw = idx_chunks[idx_chunk][idx];

          // complex frequency
          Complex w = wval + ieps;

          // build matrix and solve
          SMATRIX w_r = -w * w * in_r + ke_r;

          // if including attenuation
          if (params.attenuation()) {
            w_r += (myi + twodivpi * std::log(wval / w0)) * ke_r_atten;
          }

          // compress and find solution
          w_r.makeCompressed();
          solver.compute(w_r);
          MATRIX vec_x = solver.solve(f_r);

          // find response
          auto tmp = vec_x.block(lowidx, 0, lenidx, 1);
          auto resval = (vec_red_z.transpose() * tmp).sum();

          // add to response vector
          for (int idxr = 0; idxr < num_rec; ++idxr) {
            auto idxpl = 3 * idxr;

#pragma omp critical(torvecadd)
            {
              vec_raw(idxpl, idxw) += resval;
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
          SMATRIX H_tor = sem.H_TK(idxl).cast<Complex>();
          SMATRIX P_tor = sem.P_TK(idxl).cast<Complex>();
          SMATRIX H_tor_atten = sem.H_TA(idxl).cast<Complex>();
          H_tor.makeCompressed();
          P_tor.makeCompressed();
          H_tor_atten.makeCompressed();

          ///////////////////////////////////////////////////////////////////////
          // setting up forces and receiver vectors
          MATRIX F_VALS = sem.CalculateForce_RED_Coefficients_T(cmt, idxl, 0.0);
          MATRIX RED_C = RV_VALS * F_VALS;
          MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
          MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);

          // indices of reduced matrices at each frequency
          auto ridxtor = SpectralTools::AllIndices_TOR(
              sem, inp_model, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i2 = idx_chunks[idx_chunk].back();
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l = MATRIX::Zero(3 * params.num_receivers(),
                                          freq_chunks[idx_chunk].size());

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            auto idxw = idx_chunks[idx_chunk][idx];
            Complex w = wval + ieps;

            // lower element at this frequency
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;

            // setup
            SMATRIX mat_w_tor_red =
                H_tor.block(ridx, ridx, len_ms, len_ms) -
                w * w * P_tor.block(ridx, ridx, len_ms, len_ms);

            // if including attenuation
            if (params.attenuation()) {
              mat_w_tor_red += (myi + twodivpi * std::log(wval / w0)) *
                               H_tor_atten.block(ridx, ridx, len_ms, len_ms);
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

            MATRIX vec_sol = solver1.solve(f_red);
            auto lidx = lowidx - ridx;

            // compute responses
            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 2);

            vec_raw_l.col(idx) +=
                RED_C.cwiseProduct(testmult_red).rowwise().sum();
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

        auto t_rv_vals = clock::now();
        MATRIX RV_VALS = srinfo.RV_RED_SPH(idxl);   // receiver
        auto t_rv_vals_end = clock::now();
        tot_l_rv_vals += duration_cast<microseconds>(t_rv_vals_end - t_rv_vals);

        for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {

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

          auto t_mat = clock::now();

          // get matrices
          SMATRIX p_s = sem.P_S(idxl).cast<Complex>();
          SMATRIX h_s = sem.H_S(idxl).cast<Complex>();
          SMATRIX h_sa = sem.H_SA(idxl).cast<Complex>();
          p_s.makeCompressed();
          h_s.makeCompressed();
          h_sa.makeCompressed();
          auto t_mat_end = clock::now();
          tot_l_rmat += duration_cast<microseconds>(t_mat_end - t_mat);

          //////////////////////////////////////////////////////////////////////
          // setting up forces and receiver vectors
          auto t_f_vals = clock::now();
          MATRIX F_VALS =
              sem.CalculateForce_RED_Coefficients(cmt, idxl, 0.0);   // source
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

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            auto idxw = idx_chunks[idx_chunk][idx];

            // lower element at this frequency
            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;

            // setup
            auto t_rmat = clock::now();
            Complex w = wval + ieps;
            SMATRIX w_s = h_s.block(ridx, ridx, len_ms, len_ms) -
                          w * w * p_s.block(ridx, ridx, len_ms, len_ms);

            // if including attenuation
            if (params.attenuation()) {
              w_s += (myi + twodivpi * std::log(wval / w0)) *
                     h_sa.block(ridx, ridx, len_ms, len_ms);
            }
            w_s.makeCompressed();
            auto t_rmat_end = clock::now();
            tot_l_rmat += duration_cast<microseconds>(t_rmat_end - t_rmat);
            // std::cout << "Test f.4\n";
            auto t_f_red = clock::now();
            MATRIX f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto t_f_red_end = clock::now();
            tot_l_f_vals += duration_cast<microseconds>(t_f_red_end - t_f_red);

            auto t_factorize = clock::now();
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;

            // factorize
            if ((idxn % nskip) == 0) {
              solver1.compute(w_s);
            } else {
              solver1.factorize(w_s);
            }

            auto t_factorize_end = clock::now();
            tot_l_factorize +=
                duration_cast<microseconds>(t_factorize_end - t_factorize);

            auto t_solve = clock::now();
            MATRIX vec_sol = solver1.solve(f_red);
            auto t_solve_end = clock::now();
            tot_l_solve += duration_cast<microseconds>(t_solve_end - t_solve);
            auto lidx = lowidx - ridx;

            // compute responses
            auto t_matmult = clock::now();
            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);

            vec_raw_l.col(idx) +=
                RED_C.cwiseProduct(testmult_red).rowwise().sum();
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
  for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
    for (int idx = 0; idx < idx_chunks[idx_chunk].size(); ++idx) {
      auto wval = freq_chunks[idx_chunk][idx];
      auto idxw = idx_chunks[idx_chunk][idx];
      Complex w = wval + ieps;
      Complex mfact;
      if (params.output_type() == 0) {
        mfact = -myi / w;
      } else if (params.output_type() == 1) {
        mfact = 1.0;
      } else if (params.output_type() == 2) {
        mfact = myi * w;
      }
      vec_raw.col(idxw) *= mfact;
    }
  }
  return vec_raw;
};   // namespace SPARSESPEC
};   // namespace SPARSESPEC

#endif