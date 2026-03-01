#include <iostream>
#include <PlanetaryModel/All>
#include <new_coupling/Timer>
#include "sem_spheroidal.h"
#include "../SpectraSolver/SpectraSolver/FF"
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/filter_base.h"
#include "../SpectraSolver/SpectraSolver/src/ODE_Spectra/postprocessfunctions.h"

int
main() {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int nfreq = 5, nsolid = 1, pn = 6, lval, NQ = 5;
  bool toaug;
  double sigshift, maxstep, maxfreq = 0.002, theta_S = 0.4, phi_S = 3.1,
                            twopi = 2.0 * 3.1415926535897932;
  std::string pathpert,
      pathtoprem = "../mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";

  //////////////////////////////
  // input parameters
  std::cout << "Enter maxstep: \n";
  std::cin >> maxstep;
  std::cout << "Enter shift:\n";
  std::cin >> sigshift;
  std::cout << "Enter l:\n";
  std::cin >> lval;
  std::cout << "Enter # of frequencies:\n";
  std::cin >> nfreq;
  std::cout << "Enter solid layer:\n";
  std::cin >> nsolid;
  std::cout << "To aug or not?\n";
  std::cin >> toaug;
  std::cout << "Path to perturbed model:\n";
  std::cin >> pathpert;
  std::cout << "maxfreq:\n";
  std::cin >> maxfreq;

  //////////////////////////////
  // frequency solver parameters
  double f1 = 0.2, f2 = 1000 * maxfreq, dt = 20, tout = 100, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = 128;
  int qex = 4;
  Complex myi(0.0, 1.0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  auto prem = EarthModels::ModelInput(pathtoprem);
  SpectraSolver::FreqFull myff(f1, f2, dt, tout, df0, wtb, t1, t2, qex,
                               prem.TimeNorm());
  EarthMesh::RadialMesh _mesh(prem, NQ, 1.0, maxstep, false);
  Spheroidal::sem sem(prem, maxstep, NQ, lval + 1, 1);
  auto cmt = SourceInfo::EarthquakeCMT("./examples/bolivia_cmt_event");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // we find the spectrum in this section
  // frequencies to evaluate
  auto vec_w = myff.w();
  Complex ieps = myff.ep() * Complex(0.0, -1.0);
  MATRIX vec_raw = MATRIX::Zero(3, vec_w.size());

  // matrices and vectors for lval
  timer1.start();
  SMATRIX mat_ke = sem.MAT_KE(lval).cast<Complex>();
  SMATRIX mat_inertia = sem.MAT_IN(lval).cast<Complex>();
  mat_ke.makeCompressed();
  mat_inertia.makeCompressed();
  timer1.stop("Toroidal::sem GetStiffnessMatrix and InertiaMatrix");

  // receiver vectors
  MATRIX RV_Z = sem.RV_Z(theta_S, phi_S, lval);
  MATRIX RV_THETA = sem.RV_THETA(theta_S, phi_S, lval);
  MATRIX RV_PHI = sem.RV_PHI(theta_S, phi_S, lval);

  // force vector
  MATRIX vec_force = sem.CalculateForce(cmt, lval);

  // lu solver
  SLU solver;
  timer1.start();
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
    vec_raw(0, idx) = -w * w * RV_Z.cwiseProduct(vec_x).sum();
    vec_raw(1, idx) = -w * w * RV_THETA.cwiseProduct(vec_x).sum();
    vec_raw(2, idx) = -w * w * RV_PHI.cwiseProduct(vec_x).sum();
  };
  timer1.stop("Spheroidal::sem Solve");

  // process responses
  auto vec_r2t_b = processfunctions::filtfreq2time(vec_raw, myff);
  auto a_filt = processfunctions::fulltime2freq(vec_r2t_b, myff);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  // outputting result
  std::string pathtofile = "./work/spheroidal/groundresponse.out";
  auto file = std::ofstream(pathtofile);
  double nval = 1 / (prem.TimeNorm());
  for (int idx = 0; idx < myff.i2() + 100; ++idx) {
    file << std::setprecision(16)
         << vec_w[idx] * nval * 1000 / (2.0 * 3.1415926535) << ";"
         << a_filt(0, idx).real() << ";" << a_filt(0, idx).imag() << ";"
         << std::abs(a_filt(0, idx)) << ";" << a_filt(1, idx).real() << ";"
         << a_filt(1, idx).imag() << ";" << std::abs(a_filt(1, idx)) << ";"
         << a_filt(2, idx).real() << ";" << a_filt(2, idx).imag() << ";"
         << std::abs(a_filt(2, idx)) << std::endl;
  }
  file.close();

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
            sem.LocalToGlobal(neig, idxe, idxn);
        std::cout << "neig, idxe, idxn: " << neig << " " << idxe << " " << idxn
                  << " -> " << globidx << "\n";
      }
    }
  }
  */