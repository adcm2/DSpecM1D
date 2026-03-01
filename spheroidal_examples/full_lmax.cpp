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

template <typename FLOAT> class prem_norm {
public:
  FLOAT LengthNorm() const { return _length_norm; };
  FLOAT MassNorm() const { return _mass_norm; };
  FLOAT TimeNorm() const { return _time_norm; };

  prem_norm() {};

private:
  FLOAT _length_norm = 1000.0;
  FLOAT _mass_norm = 5515 * std::pow(_length_norm, 3.0);
  FLOAT _time_norm = 1 / std::sqrt(3.1415926535897932 * 6.67230e-11 * 5515);
};

int
main() {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read all parameters from the input file
  InputParameters params("../YSpec/input.txt");

  // earth model
  std::string cpath = params.earth_model();
  std::cout << "Current path: " << cpath << "\n";
  std::string earth_model_path = "../YSpec/" + params.earth_model();
  std::cout << "Using Earth model: " << earth_model_path << "\n";
  std::cout << "Receiver lat and lon: " << params.receivers()[0].first << ", "
            << params.receivers()[0].second << "\n";
  std::cout << "Receiver lat and lon: " << params.receivers()[1].first << ", "
            << params.receivers()[1].second << "\n";

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int nfreq = 5, nsolid = 1, lval = params.lmax(), NQ = 5;
  bool toaug = false;
  double maxstep = 0.05;
  const double twopi = 2.0 * 3.1415926535897932;
  std::string pathpert;

  //////////////////////////////
  // frequency solver parameters
  double f1 = params.fmin_mHz(), f2 = params.fmax_mHz(),
         dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout + 1.0;
  int qex = 4;
  Complex myi(0.0, 1.0);
  // lval = params.lmax();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  // auto prem = EarthModels::ModelInput(pathtoprem);

  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, true);
  // std::cout << "Earth model loaded.\n";
  // std::cout << "Length norm: " << prem.LengthNorm() << " m\n";
  // std::cout << "Mass norm:   " << prem.MassNorm() << " kg\n";
  // std::cout << "Time norm:   " << prem.TimeNorm() << " s\n\n";

  // frequency class
  SpectraSolver::FreqFull myff(params.f11(), params.f12(), params.f21(),
                               params.f22(), dt, tout, df0, wtb, t1, t2, qex,
                               prem.TimeNorm());

  // std::cout << "Frequency setup:\n";
  // mesh
  // EarthMesh::RadialMesh _mesh(prem, NQ, 1.0, maxstep, false);

  // sem class
  Full1D::sem sem(prem, maxstep, NQ, lval);

  // std::cout << "SEM model created successfully." << std::endl;

  // source information
  auto cmt = SourceInfo::EarthquakeCMT(params);
  // std::cout << "cmt read in\n";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // compare sem with transverse isotropy and purely isotropic
  //  Spheroidal::sem sem2(prem, maxstep, NQ, lval, 1);

  // compare sem ke matrices
  // SMATRIX mat_ke1 = sem.MAT_KE(1).cast<Complex>();
  // SMATRIX mat_ke2 = sem2.MAT_KE(1).cast<Complex>();
  // go through matrices and compare the values
  // for (int k = 0; k < mat_ke1.outerSize(); ++k) {
  //   for (typename SMATRIX::InnerIterator it1(mat_ke1, k), it2(mat_ke2, k);
  //        it1 && it2; ++it1, ++it2) {
  //     if (it1.row() != it2.row() || it1.col() != it2.col()) {
  //       std::cerr << "Error: matrices have different structure at ("
  //                 << it1.row() << "," << it1.col() << ")\n";
  //       return 1;
  //     }
  //     if (std::abs((it1.value() - it2.value()) /
  //                  std::max(std::abs(it1.value()), std::abs(it2.value()))) >
  //         1e-10) {
  //       std::cerr << "Error: matrices differ at (" << it1.row() << ","
  //                 << it1.col() << "): " << it1.value() << " vs " <<
  //                 it2.value()
  //                 << "\n";
  //       return 1;
  //     }
  //   }
  // }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex ieps = myff.ep() * Complex(0.0, -1.0);
  // std::cout << "\neps: " << myff.ep() << "\n\n";
  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  timer1.start();
  SLU solver, solverT;

  auto idxout = sem.LtG_S(0, 25, NQ - 1);
  double momentnorm =
      prem.MassNorm() * std::pow(prem.LengthNorm() / prem.TimeNorm(), 2.0);
  // std::cout << "Moment norm: " << momentnorm << "\n\n";
  std::vector<MATRIX> vec_response_save(lval);

  // std::cout << params.lmin() << " to " << lval << "\n";

  ///////////////////////////////////
  // getting minimum and maximum l values
  bool radial_mode = (params.lmin() == 0);
  int lmin_s = std::max(params.lmin(), 1);
  int lmax_s = params.lmax();
  int lmin_t = std::max(params.lmin(), 2);
  int lmax_t = lmax_s;
  ///////////////////////////////////

  timer1.start();
  // do spheroidals:

  for (int idxl = lmin_s; idxl < lmax_s + 1; ++idxl) {
    SMATRIX mat_ke = sem.MAT_KE(idxl).cast<Complex>();
    SMATRIX mat_inertia = sem.MAT_IN(idxl).cast<Complex>();
    mat_ke.makeCompressed();
    mat_inertia.makeCompressed();

    // analyze pattern for quicker lu decomposition
    {
      SMATRIX mat_test = mat_ke + mat_inertia;
      solver.analyzePattern(mat_test);
    }

    // calculate force vector
    MATRIX vec_force = sem.CalculateForce(cmt, idxl);
    // std::cout << "Force values: \n"
    //           << vec_force.block(idxout - 1, 0, 2, 3) << "\n\n";

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
      solver.factorize(mat_w);
      MATRIX vec_x = solver.solve(vec_fw);
      if (idx * myff.df() < 0.232 / 1000 * prem.TimeNorm()) {
        // std::cout << idx * myff.df() / prem.TimeNorm() * 1000 << "\n";
        vec_response_save[idxl - 1] = vec_x;
      }
      // std::cout << "\n";
      // compute responses
      for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {

        // receiver vectors
        auto RV_Z = sem.RV_Z(params, idxl, idxr);
        auto RV_THETA = sem.RV_THETA(params, idxl, idxr);
        auto RV_PHI = sem.RV_PHI(params, idxl, idxr);

        // index
        auto idxpl = 3 * idxr;

        // find response
        vec_raw(idxpl, idx) -= w * w * RV_Z.cwiseProduct(vec_x).sum();
        vec_raw(idxpl + 1, idx) += w * w * RV_THETA.cwiseProduct(vec_x).sum();
        vec_raw(idxpl + 2, idx) -= w * w * RV_PHI.cwiseProduct(vec_x).sum();

        // NOTE: signs are according to the convention used in YSpec. In
        // particular for theta its reversed so it is the north componenet
      }
    };
  };
  timer1.stop("Spheroidal Solve");

  timer1.start();
  /*
    // do toroidals:
    for (int idxl = lmin_t; idxl < lmax_t + 1; ++idxl) {
      SMATRIX mat_ke = sem.MAT_KE_T(idxl).cast<Complex>();
      SMATRIX mat_inertia = sem.MAT_IN_T(idxl).cast<Complex>();
      mat_ke.makeCompressed();
      mat_inertia.makeCompressed();
      // calculate force vector
      MATRIX vec_force = sem.CalculateForce_T(cmt, idxl);
      // std::cout << "Force values: \n"
      //           << vec_force.block(idxout - 1, 0, 2, 3) << "\n\n";
      // std::cout << "Rows in ke: " << mat_ke.rows()
      //           << " Cols in ke: " << mat_ke.cols() << "\n";
      // std::cout << "Rows in inertia: " << mat_inertia.rows()
      //           << " Cols in inertia: " << mat_inertia.cols() << "\n";
      // std::cout << "Size of force vector: " << vec_force.rows() << " x "
      //           << vec_force.cols() << "\n";
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
          auto RV_THETA = sem.RV_THETA_T(params, idxl, idxr);
          auto RV_PHI = sem.RV_PHI_T(params, idxl, idxr);

          // index
          auto idxpl = 3 * idxr;

          // find response
          vec_raw(idxpl + 1, idx) += w * w * RV_THETA.cwiseProduct(vec_x).sum();
          vec_raw(idxpl + 2, idx) -= w * w * RV_PHI.cwiseProduct(vec_x).sum();

          // NOTE: signs are according to the convention used in YSpec. In
          // particular for theta its reversed so it is the north componenet
        }
      };
    }*/

  timer1.stop("Toroidal Solve");

  // normalise
  auto accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  vec_raw *= accel_norm;

  // process responses
  auto vec_r2t_b = processfunctions::filtfreq2time(vec_raw, myff);
  auto a_filt = processfunctions::fulltime2freq(vec_r2t_b, myff);
  // auto vec_filt_t = processfunctions::filtfreq2time(a_filt, myff, false);
  // auto a_filt_final = processfunctions::fulltime2freq(vec_filt_t, myff);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // testing read in of yspec
  // std::cout << "First output\n";
  std::string yspec_path = "../YSpec/output/yspec.out.1";
  YSPECREADER::DataColumns yspec_data(yspec_path);

  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (int idx = 0; idx < yspec_data.getColumn1().size(); ++idx) {
    yspec_t(0, idx) += yspec_data.getColumn2()[idx];
    yspec_t(1, idx) += yspec_data.getColumn3()[idx];
    yspec_t(2, idx) += yspec_data.getColumn4()[idx];
  }
  auto a_filt_yspec = processfunctions::fulltime2freq(yspec_t, myff);

  // read in mineos output
  std::string mineos_path = "../mineos-1.0.2/DEMO/MYEX/Syndat_ASC_NOHEADER/"
                            "Syndat.2000014:23:37:10.TLY.LHZ.ASC";
  MINEOSREADER::DataColumns mineos_data(mineos_path);
  std::string mineos_path1 = "../mineos-1.0.2/DEMO/MYEX/Syndat_ASC_NOHEADER/"
                             "Syndat.2000014:23:37:10.TLY.LHN.ASC";
  MINEOSREADER::DataColumns mineos_data1(mineos_path1);
  std::string mineos_path2 = "../mineos-1.0.2/DEMO/MYEX/Syndat_ASC_NOHEADER/"
                             "Syndat.2000014:23:37:10.TLY.LHE.ASC";
  MINEOSREADER::DataColumns mineos_data2(mineos_path2);

  Eigen::MatrixXd mineos_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (int idx = 0; idx < mineos_data.getColumn1().size(); ++idx) {
    mineos_t(0, idx) += mineos_data.getColumn2()[idx] * 1e-9;
    mineos_t(1, idx) += mineos_data1.getColumn2()[idx] * 1e-9;
    mineos_t(2, idx) += mineos_data2.getColumn2()[idx] * 1e-9;
  }
  auto a_mineos_0 = processfunctions::fulltime2freq(mineos_t, myff);
  auto vec_filt_t_mineos =
      processfunctions::filtfreq2time(a_mineos_0, myff, false);
  auto a_filt_mineos = processfunctions::fulltime2freq(vec_filt_t_mineos, myff);
  // auto a_filt_mineos =

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  // outputting result
  // std::cout << "Second output\n";
  std::string pathtofile = "./work/spheroidal/groundresponse_l_transverse.out";
  std::ofstream file(pathtofile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  const auto &vec_w_ref = vec_w;             // alias
  const std::size_t nw = vec_w_ref.size();   // use actual vector size

  // Write header (optional)
  // file <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file.setf(std::ios::fixed);
  file << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has compatible
    // columns
    file << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
         << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
         << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
         << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
         << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
         << std::abs(a_filt(2, idx)) << ";" << a_filt_yspec(0, idx).real()
         << ";" << a_filt_yspec(0, idx).imag() << ";"
         << std::abs(a_filt_yspec(0, idx)) << ";" << a_filt_yspec(1, idx).real()
         << ";" << a_filt_yspec(1, idx).imag() << ";"
         << std::abs(a_filt_yspec(1, idx)) << ";" << a_filt_yspec(2, idx).real()
         << ";" << a_filt_yspec(2, idx).imag() << ";"
         << std::abs(a_filt_yspec(2, idx)) << ";"
         << a_filt_mineos(0, idx).real() << ";" << a_filt_mineos(0, idx).imag()
         << ";" << std::abs(a_filt_mineos(0, idx)) << ";"
         << a_filt_mineos(1, idx).real() << ";" << a_filt_mineos(1, idx).imag()
         << ";" << std::abs(a_filt_mineos(1, idx)) << ";"
         << a_filt_mineos(2, idx).real() << ";" << a_filt_mineos(2, idx).imag()
         << ";" << std::abs(a_filt_mineos(2, idx)) << '\n';
  }
  file.close();

  //////////////////////////////////////////////////////////////////////////////
  // Repeat output for raw data
  // std::cout << "Third output\n";
  std::string pathtofile2 =
      "./work/spheroidal/groundresponse_raw_transverse.out";
  std::ofstream file2(pathtofile2);
  if (!file2) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  // Write header (optional)
  // file <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file2.setf(std::ios::fixed);
  file2 << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has compatible
    // columns
    file2 << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
          << vec_raw(0, idx).real() << ';' << vec_raw(0, idx).imag() << ';'
          << std::abs(vec_raw(0, idx)) << ';' << vec_raw(1, idx).real() << ';'
          << vec_raw(1, idx).imag() << ';' << std::abs(vec_raw(1, idx)) << ';'
          << vec_raw(2, idx).real() << ';' << vec_raw(2, idx).imag() << ';'
          << std::abs(vec_raw(2, idx)) << '\n';
  }
  file2.close();

  //////////////////////////////////////////////////////////////////////////////
  // Repeat output for raw data
  // std::cout << "Fourth output\n";
  std::string pathtofile3 =
      "./work/spheroidal/groundresponse_time_transverse.out";
  std::ofstream file3(pathtofile3);
  if (!file3) {
    std::cerr << "Error: unable to open output file: " << pathtofile3 << "\n";
    return 1;
  }

  // Write header (optional)
  // file <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file3.setf(std::ios::fixed);
  file3 << std::setprecision(16);

  for (std::size_t idx = 0; idx < vec_r2t_b.cols(); ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has compatible
    // columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file3 << (idx * myff.dt() * prem.TimeNorm()) << ';' << vec_r2t_b(0, idx)
          << ';' << vec_r2t_b(1, idx) << ';' << vec_r2t_b(2, idx) << '\n';
    if (tval > params.t_out() * 60.0) {
      break;
    }
  }
  file3.close();

  //////////////////////////////////////////////////////////////////////////////
  // Repeat output for raw data
  // std::cout << "Fifth output\n";
  std::string pathtofile4 =
      "./work/spheroidal/groundresponse_t_yspec_transverse.out";
  std::ofstream file4(pathtofile4);
  if (!file4) {
    std::cerr << "Error: unable to open output file: " << pathtofile4 << "\n";
    return 1;
  }
  file4.setf(std::ios::fixed);
  file4 << std::setprecision(16);

  for (std::size_t idx = 0; idx < yspec_data.getColumn1().size(); ++idx) {
    file4 << yspec_data.getColumn1()[idx] << ";" << yspec_data.getColumn2()[idx]
          << ";" << yspec_data.getColumn3()[idx] << ";"
          << yspec_data.getColumn4()[idx] << "\n";
  }
  file4.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string pathtofile_mineos =
      "./work/spheroidal/groundresponse_t_mineos_transverse.out";
  std::ofstream file_mineos(pathtofile_mineos);
  if (!file_mineos) {
    std::cerr << "Error: unable to open output file: " << pathtofile_mineos
              << "\n";
    return 1;
  }
  file_mineos.setf(std::ios::fixed);
  file_mineos << std::setprecision(16);
  for (std::size_t idx = 0; idx < vec_filt_t_mineos.cols(); ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has compatible
    // columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_mineos << (idx * myff.dt() * prem.TimeNorm()) << ';'
                << vec_filt_t_mineos(0, idx) << ';' << vec_filt_t_mineos(1, idx)
                << ';' << vec_filt_t_mineos(2, idx) << ";" << vec_r2t_b(0, idx)
                << ";" << vec_r2t_b(1, idx) << ";" << vec_r2t_b(2, idx) << '\n';
    if (tval > params.t_out() * 60.0) {
      break;
    }
  }
  file_mineos.close();

  //////////////////////////////////////////////////////////////////////////////
  // Repeat output for raw data
  // std::cout << "Sixth output\n";
  // std::string pathtofile5 = "./work/spheroidal/fullresponse_single_w.out";
  // std::ofstream file5(pathtofile5);
  // if (!5) {
  //   std::cerr << "Error: unable to open output file: " << pathtofile5 <<
  //   "\n"; return 1;
  // }
  // file5.setf(std::ios::fixed);
  // file5 << std::setprecision(16);

  // for (auto idxe = 0; idxe < sem.mesh().NE(); ++idxe) {
  //   for (auto idxq = 0; idxq < NQ; ++idxq) {
  //     auto uidx = sem.LtG_S(0, idxe, idxq);
  //     auto vidx = sem.LtG_S(1, idxe, idxq);
  //     auto pidx = sem.LtG_S(2, idxe, idxq);
  //     auto rval = sem.mesh().NodeRadius(idxe, idxq);
  //     auto uval = vec_response_save[0](uidx, 0);
  //     auto vval = vec_response_save[0](vidx, 0);
  //     auto pval = vec_response_save[0](pidx, 0);
  //     file5 << rval << ";" << uval.real() << ";" << uval.imag() << ";"
  //           << vval.real() << ";" << vval.imag() << ";" << pval.real() << ";"
  //           << pval.imag() << "\n";
  //   }
  // }
  // file5.close();

  return 0;
}

/*
  // find fs boundaries
  auto fsb = _mesh.FS_Boundaries();
  int idxnum = 1;
  for (auto v : fsb) {
    std::cout << "FS boundary #" << idxnum++ << ":\n"
              << "Lower element: " << v
              << " with upper radius: " << _mesh.EUR(v) << " and is solid? "
              << _mesh.IsSolid(v) << "\n"
              << "Upper element: " << v + 1
              << " with lower radius: " << _mesh.ELR(v + 1) << " and is solid? "
              << _mesh.IsSolid(v + 1) << "\n\n";
  }

  std::vector<int> vec_offset{-2};
  {
    std::size_t totnum = 0;
    for (int idx = 0; idx < _mesh.NE(); ++idx) {
      auto tmp = vec_offset[idx];
      if (idx == fsb[totnum]) {
        tmp += 1;
        totnum += 1;
      }
      vec_offset.push_back(tmp);
    }
  }
  int idxoutput = 0;
  for (auto v : vec_offset) {
    std::cout << idxoutput++ << " " << v << "\n ";
  }
  std::cout << "\n";
  auto local_to_global = [&fsb, &vec_offset, NQ](int neig, int idx_e,
                                                 int idx_n) {
    std::size_t retval;
    if ((idx_e == 0) && (idx_n == 0)) {
      if (neig == 2) {
        retval = 0;
      } else {
        std::cout << "Error: idx_e = 0 but neig != 2\n";
      }
    } else {
      int offset_val = vec_offset[idx_e];
      if (idx_n == 0) {
        if ((std::find(fsb.begin(), fsb.end(), idx_e - 1) != fsb.end())) {
          if (neig == 1) {
            offset_val += 1;
          } else {
            offset_val -= 1;
          }
        }
      }
      retval = (3 * idx_e * (NQ - 1) + idx_n * 3 + neig) + offset_val;
    }

    return retval;
  };
  std::cout << "Check local to global:\n";
  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    for (int idxn = 0; idxn < NQ; ++idxn) {
      for (int idxg = 0; idxg < 3; ++idxg) {
        std::cout << "neig, idxe, idxn: " << idxg << " " << idxe << " " << idxn
                  << " -> " << local_to_global(idxg, idxe, idxn) << "\n";
      }
    }
  }
  std::cout << "\n";

  std::cout << "Spheroidal SEM model created successfully." << std::endl;
  std::cout << "Testing local to global mapping:" << std::endl;
  for (int idxe = 0; idxe < sem.mesh().NE(); ++idxe) {
    for (int idxn = 0; idxn < NQ; ++idxn) {
      for (int neig = 0; neig < 3; ++neig) {
        std::size_t globidx =
            sem.LtG_S(neig, idxe, idxn);
        std::cout << "neig, idxe, idxn: " << neig << " " << idxe << " " << idxn
                  << " -> " << globidx << "\n";
      }
    }
  }
  */