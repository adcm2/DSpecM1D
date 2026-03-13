#ifndef FULL_SPEC_GUARD_H
#define FULL_SPEC_GUARD_H

#include <iostream>
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include "ReadStation.h"
#include "SourceInfo.h"
#include "StartElement.h"
#include "SpectraMaster.h"
#include "ParamInfo.h"
#include <omp.h>
#include "FEM_Preconditioner.h"
#include "BiCGSTABT.h"
#include "ParamRedInfo.h"
#include "SR_Info.h"
#include "SpecHelpers.h"

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
  SpecConstants sc(myff.ep(), inp_model.TREF());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
  auto twodivpi = sc.twodivpi;

  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  SLU solver, solver1;

  ///////////////////////////////////
  // getting minimum and maximum l values
  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamInfo param_info(params, lmax);
  auto NQ = sem.mesh().NN();

  // decide which modes to include
  auto mtype = params.type();
  ModeFlags flags = resolveModeFlags(mtype, lmin, lmax);
  bool inc_rad = flags.inc_rad;
  bool inc_tor = flags.inc_tor;
  bool inc_sph = flags.inc_sph;

  // change lmin for toroidal and spheroidal modes
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();

  ///////////////////////////////////
  auto idx_source = sem.Source_Element(cmt);

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
          w_r += attenFactor(vec_w[idx], w0, twodivpi, myi) * ke_r_atten;
        }
        w_r.makeCompressed();
        solver.compute(w_r);
        MATRIX vec_x = solver.solve(f_r);

        // compute responses
        for (int idxr = 0; idxr < num_rec; ++idxr) {
          auto idxpl = 3 * idxr;
          Complex mfact = 1.0;
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
        auto ridxtor =
            SpectralTools::AllIndices_TOR(sem, idxl, myff, idx_source, nskip);
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
            mat_w_tor_red += attenFactor(vec_w[idx], w0, twodivpi, myi) *
                             H_tor_atten.block(ridx, ridx, len_ms, len_ms);
          }
          mat_w_tor_red.makeCompressed();

          auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
          auto idxn = myff.i2() - idx - 1;
          factorizeOrCompute(solver1, mat_w_tor_red, idxn, nskip);

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

        auto vec_ridx =
            SpectralTools::AllIndices_SPH(sem, idxl, myff, idx_source, nskip);

        ///////////////////////////////////////////////////////////////////////////////
        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          // indices
          std::size_t idx_rs = vec_ridx[idx - myff.i1()];
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
          factorizeOrCompute(solver1, mat_sph, idxn, nskip);

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
    Complex w = vec_w[idx] + ieps;
    Complex mfact = outputFactor(params.output_type(), w, myi);
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
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  auto vec_w = myff.w();
  SpecConstants sc(myff.ep(), inp_model.TREF());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
  auto twodivpi = sc.twodivpi;
  auto twopid = sc.twopid;
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

  auto mtype = params.type();
  ModeFlags flags = resolveModeFlags(mtype, lmin, lmax);
  bool inc_rad = flags.inc_rad;
  bool inc_tor = flags.inc_tor;
  bool inc_sph = flags.inc_sph;

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
            w_r += attenFactor(wval, w0, twodivpi, myi) * ke_r_atten;
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
              sem, idxl, freq_chunks[idx_chunk], idx_source, nskip);
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
              mat_w_tor_red += attenFactor(wval, w0, twodivpi, myi) *
                               H_tor_atten.block(ridx, ridx, len_ms, len_ms);
            }
            mat_w_tor_red.makeCompressed();
            auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;
            factorizeOrCompute(solver1, mat_w_tor_red, idxn, nskip);

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

#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {

        MATRIX RV_VALS = srinfo.RV_RED_SPH(idxl);

        for (auto idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {

          Full1D::specsem &sem = sems[idx_chunk];
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
          auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;

          // get matrices
          SMATRIX p_s = sem.P_S(idxl).cast<Complex>();
          SMATRIX h_s = sem.H_S(idxl).cast<Complex>();
          SMATRIX h_sa = sem.H_SA(idxl).cast<Complex>();
          p_s.makeCompressed();
          h_s.makeCompressed();
          h_sa.makeCompressed();

          // setting up forces and receiver vectors
          MATRIX F_VALS = sem.CalculateForce_RED_Coefficients(cmt, idxl, 0.0);
          MATRIX RED_C = RV_VALS * F_VALS;
          MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
          MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);

          auto ridxsph = SpectralTools::AllIndices_SPH(
              sem, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i2 = idx_chunks[idx_chunk].back();
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l = MATRIX::Zero(3 * params.num_receivers(),
                                          freq_chunks[idx_chunk].size());

          for (int idx = idx_chunks[idx_chunk].size() - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];

            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;

            Complex w = wval + ieps;
            SMATRIX w_s = h_s.block(ridx, ridx, len_ms, len_ms) -
                          w * w * p_s.block(ridx, ridx, len_ms, len_ms);

            if (params.attenuation()) {
              w_s += attenFactor(wval, w0, twodivpi, myi) *
                     h_sa.block(ridx, ridx, len_ms, len_ms);
            }
            w_s.makeCompressed();

            MATRIX f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            auto idxn = idx_chunks[idx_chunk].size() - idx - 1;
            factorizeOrCompute(solver1, w_s, idxn, nskip);

            MATRIX vec_sol = solver1.solve(f_red);
            auto lidx = lowidx - ridx;

            auto testmult_red = RV_BASE * vec_sol.block(lidx, 0, lenidx, 4);
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
      Complex mfact = outputFactor(params.output_type(), w, myi);

      vec_raw.col(idxw) *= mfact;
    }
  }
  return vec_raw;
};   // namespace SPARSESPEC
};   // namespace SPARSESPEC

#endif