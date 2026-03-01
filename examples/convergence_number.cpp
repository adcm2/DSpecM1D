#include <Eigen/Core>
#include <Eigen/LU>
#include <GaussQuad/All>
#include <new_coupling/All>
#include <PlanetaryModel/All>

#include "convergence_number.h"

int
main() {
  std::cout << "Hello\n";
  Timer timer1;
  timer1.start();

  // paths
  std::string pathtopremeig =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_D";
  std::string pathtoprempert =
      "/space/adcm2/mineos-1.0.2/OUTPUT/eprem_noocean_S_D_perturb";
  // std::string pathtoprem =
  //     "/space/adcm2/mineos-1.0.2/DEMO/models/prem_noocean_noattenuation.txt";
  // std::string pathtopremperturb = "/space/adcm2/mineos-1.0.2/DEMO/models/"
  //                                 "prem_noocean_noattenuation_perturb.txt";

  std::string pathtoprem =
      "/home/adcm2/space/mineos-1.0.2/DEMO/models/prem.200.no.noatten.txt";
  std::string pathtopremperturb = "/home/adcm2/space/mineos-1.0.2/DEMO/models/"
                                  "prem.200.no.noatten.perturb.txt";

  // models
  // auto prem = TestTools::EarthModel(pathtoprem);
  // auto premperturb = TestTools::EarthModel(pathtopremperturb);
  auto prem = EarthModels::ModelInput(pathtoprem);
  auto premperturb = EarthModels::ModelInput(pathtopremperturb);

  // std::cout << "# elements: " << testradmesh.NumberOfElements()
  //           << ", outer: " << testradmesh.OuterRadius()
  //           << ", planet: " << testradmesh.PlanetRadius() << "\n";
  // std::cout << "Outer Layer: "
  //           << testradmesh.LayerNumber(testradmesh.NumberOfElements() - 1)
  //           << "\n";

  // mode catalogues
  ModeCoupling::modecataloguecontinuous modes_prem(pathtopremeig, prem);
  ModeCoupling::modecataloguecontinuous modes_perturbed(pathtoprempert,
                                                        premperturb);

  timer1.stop("Time to read in modes");
  // for (int idx = 0; idx < 5; ++idx) {
  //    auto testmode = modes_prem.singlemode(idx);
  //    std::cout << testmode.n() << " " << testmode.l() << " " << testmode.q()
  //              << " " << testmode.gv() << " " << testmode.w() << "\n";
  // }
  std::cout << "Number of modes: " << modes_prem.NumberOfModes() << "\n";

  // layers
  int numlayers = prem.NumberOfLayers();

  int maxmodes = modes_prem.NumberOfModes();
  maxmodes = 150;
  // timer1.start();
  // Eigen::MatrixXd mat_norms = mode_gram_matrix(modes_prem, prem, maxmodes);
  // timer1.stop("Time for matrix assembly");

  // testing mode summary class
  ModeCoupling::modesummary prem_summary(modes_prem, prem, 8);
  ModeCoupling::modesummary2 prem_summary2(modes_prem, prem, 10);

  std::cout << "Frequency of 199: " << prem_summary.w(0) << " "
            << prem_summary.w(199) << "\n\n\n\n";

  // for (int idx = 0; idx < prem_summary2.fullmesh().NumberOfElements();
  // ++idx) {
  //    for (int idxq = 0; idxq < 7; ++idxq) {
  //       std::cout << idx << " " << prem_summary2.fullmesh().LayerNumber(idx)
  //                 << " " << idxq << " "
  //                 << prem_summary2.fullmesh().NodeRadius(idx, idxq) << "\n";
  //    }
  // }

  timer1.start();
  Eigen::MatrixXd mat_norms = mode_gram_matrix_s(prem_summary, prem, maxmodes);
  timer1.stop("Time matrix assembly 1");
  Eigen::MatrixXd mat_norms2 =
      mode_gram_matrix_s(prem_summary, prem, maxmodes, true);

  timer1.start();
  Eigen::MatrixXd mat_norms3 =
      mode_gram_matrix_s2(prem_summary2, prem, maxmodes, true);
  timer1.stop("Time matrix assembly 2");
  // std::cout << mat_norms.block(0, 0, 10, 10) -
  //                  mat_norms_check.block(0, 0, 10, 10)
  //           << "\n";
  // std::cout << mat_norms2.block(0, 0, 5, 5) << "\n";
  // std::cout << mat_norms2.block(0, 0, 5, 5) - mat_norms3.block(0, 0, 5, 5)
  //           << "\n";
  // std::cout << mat_norms2.block(0, 0, 5, 5) - mat_norms.block(0, 0, 5, 5)
  //           << "\n";
  // std::cout << "Size: " << mat_norms2.rows() << "\n";

  timer1.start();
  auto l21_pert = modes_perturbed.singlemode(0);
  auto l21_prem = modes_prem.singlemode(0);
  timer1.stop("Time for getting modes");

  // test radial mesh
  int npointsrad = 10;
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npointsrad);
  // auto testradmesh = Radial_Tools::RadialMesh(prem, q, 0.01, 1.0, false);
  auto testradmesh = EarthMesh::RadialMesh(prem, npointsrad, 0.01, 1.0, false);

  // std::cout << "Planet: " << testradmesh.PlanetRadius()
  //           << ", outer: " << testradmesh.OuterRadius() << "\n";
  auto result1 = modescalarproduct(l21_prem, l21_prem, prem);
  auto result2 =
      modescalarproduct_2(l21_prem, l21_prem, prem, testradmesh, npointsrad);
  // std::cout << std::setprecision(16) << mat_norms2(0, 0) << " "
  //           << mat_norms3(0, 0) << " " << result1 << " " << result2 << "\n";
  // std::cout << std::setprecision(16) << mat_norms2(0, 0) << " "
  //           << mat_norms3(0, 0) << " " << result1 * l21_prem.w() *
  //           l21_prem.w()
  //           << " "
  //           << result2 * l21_prem.w() * l21_prem.w() /
  //                  (prem.OuterRadius() * prem.OuterRadius() *
  //                   prem.OuterRadius())
  //           << "\n";
  //  for (int idxe = 0; idxe < 30; ++idxe){
  //    std::cout << idxe << " " << testradmesh.ElementUpperRadius(idxe) <<" "
  //    << testradmesh.LayerNumber(idxe) << "\n";
  //  }
  // test mode force with lagrange
  // std::vector<int> lagint{0, 2};
  // double testforce =
  //     mode_force_lagrange_vector(l21_pert, prem, maxmodes, lagint);

  timer1.start();
  auto vec_force = mode_force_vector_s(prem_summary, l21_pert, prem, maxmodes);
  auto vec_force2 =
      mode_force_vector_s(prem_summary, l21_pert, prem, maxmodes, true);
  auto vec_force3 =
      mode_force_vector_s2(prem_summary2, l21_pert, prem, maxmodes, true);
  timer1.stop("Time for force vector assembly 2");

  Eigen::FullPivLU<Eigen::MatrixXd> solverf(mat_norms);
  Eigen::FullPivLU<Eigen::MatrixXd> solverf2(mat_norms2);
  Eigen::FullPivLU<Eigen::MatrixXd> solverf3(mat_norms3);
  Eigen::VectorXd vec_fullsol = solverf.solve(vec_force);
  Eigen::VectorXd vec_fullsol2 = solverf2.solve(vec_force2);
  Eigen::VectorXd vec_fullsol3 = solverf3.solve(vec_force3);

  // for (int idx = 0; idx < 10; ++idx) {
  //    std::cout << (vec_fullsol2(idx) - vec_fullsol3(idx)) /
  //    vec_fullsol2(idx)
  //              << "\n";
  // }
  // vec_fullsol2
  // std::cout << mat_norms2.block(mat_norms2.rows() - 15, mat_norms2.rows() -
  // 15,
  //                               14, 14)
  //           << "\n\n";
  // std::cout << vec_force2.tail(14) << "\n\n";

  // std::cout << vec_fullsol2.tail(14) << "\n\n";

  timer1.start();

  // Eigen::VectorXd vec_sol2 = solver.solve(vec_force2);
  std::vector<double> vec_error(maxmodes, 0.0), vec_error2(maxmodes, 0.0),
      vec_error_v(maxmodes, 0.0), vec_errorh10(maxmodes, 0.0),
      vec_errorh12(maxmodes, 0.0), vec_traction_r(maxmodes, 0.0),
      vec_traction_r_pert(maxmodes, 0.0),
      vec_traction_r_pertexact(maxmodes, 0.0);
  double normofmode = l2norm(l21_pert, premperturb);
  double normofmode_v = l2norm_v(l21_pert, premperturb);
  double h10normofmode = h10norm(l21_pert, premperturb);

  for (int idx = 1; idx < maxmodes + 1; ++idx) {
    Eigen::MatrixXd mat_small = mat_norms.block(0, 0, idx, idx);
    Eigen::FullPivLU<Eigen::MatrixXd> solver(mat_small);
    Eigen::VectorXd vec_small = vec_force.head(idx);
    Eigen::VectorXd vec_sol = solver.solve(vec_small);

    // radial mesh
    Eigen::MatrixXd mat_small2 = mat_norms3.block(0, 0, idx, idx);
    Eigen::FullPivLU<Eigen::MatrixXd> solver2(mat_small2);
    Eigen::VectorXd vec_small2 = vec_force3.head(idx);
    Eigen::VectorXd vec_sol2 = solver2.solve(vec_small2);

    // correct for frequency norm factor
    for (int idxc = 0; idxc < vec_sol.rows(); ++idxc) {
      // vec_sol(idxc) *= modes_prem.singlemodep(idxc).w() /
      //                  modes_perturbed.singlemodep(0).w();
    }
    // vec_error[idx - 1] =
    //     l2normdifference_s(vec_sol, prem_summary, l21_pert, prem) /
    //     normofmode;
    vec_error[idx - 1] =
        l2normdifference_s_full(vec_sol, prem_summary, l21_pert, prem) /
        normofmode;
    vec_error2[idx - 1] =
        l2normdifference_s_full2(vec_sol2, prem_summary2, l21_pert, prem) /
        normofmode;
    vec_error_v[idx - 1] =
        l2normdifference_sv(vec_sol, prem_summary, l21_pert, prem) /
        normofmode_v;
    vec_errorh10[idx - 1] =
        h10normdifference_s(vec_sol, prem_summary, l21_pert, prem) /
        h10normofmode;
    vec_errorh12[idx - 1] =
        h10normdifference_s2(vec_sol2, prem_summary2, l21_pert, prem) /
        h10normofmode;
    vec_traction_r[idx - 1] = traction_basis_sum(vec_sol, prem_summary, prem);
    vec_traction_r_pert[idx - 1] =
        traction_basis_sum(vec_sol, prem_summary, premperturb);
    vec_traction_r_pertexact[idx - 1] =
        traction_singlemode(l21_pert, premperturb);
  }
  timer1.stop("Time for solving at different number of basis functions");

  for (int idx = 0; idx < maxmodes - 1; ++idx) {
    if (vec_error[idx + 1] > vec_error[idx]) {
      std::cout << "Error got bigger: " << idx << "\n";
    }
  }
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  // outputting result
  std::string pathtofile = "./work/normerror.out";
  auto file = std::ofstream(pathtofile);
  for (int i = 0; i < maxmodes; ++i) {
    file << std::setprecision(16) << i + 1 << ";" << vec_error[i] << ";"
         << vec_error_v[i] << ";" << vec_errorh10[i] << ";" << vec_error2[i]
         << ";" << vec_errorh12[i] << std::endl;
  };
  file.close();

  // outputting result
  std::string pathtofile2 = "./work/prempert.out";
  auto file2 = std::ofstream(pathtofile2);
  {
    auto qx = prem_summary.q().Points();
    for (int idxlayer = 0; idxlayer < prem.NumberOfLayers(); ++idxlayer) {
      auto vec_radii = prem.LayerRadii(idxlayer);

      for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];
        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        // double tmpsum = 0.0;

        for (int idxq = 0; idxq < qx.size(); ++idxq) {
          double xscale = delta * qx[idxq] + xpl;
          file2 << std::setprecision(16) << xscale << ";"
                << l21_pert.u(idxlayer)(xscale) << ";"
                << modesum_s(vec_fullsol, prem_summary, prem, idxlayer, idxpl,
                             idxq)
                << ";" << l21_prem.u(idxlayer)(xscale) << ";"
                << modesum_s(vec_fullsol2, prem_summary, prem, idxlayer, idxpl,
                             idxq, true)
                << ";" << prem_summary.uvec(maxmodes - 1, idxlayer, idxpl, idxq)
                << std::endl;
          // std::cout << modesum_s(vec_fullsol2, prem_summary, prem,
          //                        idxlayer, idxpl, idxq, true)
          //           << "\n";
        }
      };
    };
  };

  file2.close();

  // outputting result
  std::string pathtofile10 = "./work/prempert2.out";
  auto file10 = std::ofstream(pathtofile10);
  {
    auto qx = prem_summary.q().Points();

    for (int idxe = 0; idxe < prem_summary2.fullmesh().NE(); ++idxe) {
      int idxlayer = prem_summary2.fullmesh().LayerNumber(idxe);
      for (int idxq = 0; idxq < prem_summary2.npoints(); ++idxq) {
        double x = prem_summary2.fullmesh().NodeRadius(idxe, idxq);
        file10 << std::setprecision(16) << x << ";" << l21_prem.u(idxlayer)(x)
               << ";" << l21_pert.u(idxlayer)(x) << ";"
               << modesum_s2(vec_fullsol3, prem_summary2, prem, idxe, idxq,
                             true)
               << std::endl;
      }
      // };
    };
  };

  file10.close();

  // outputting result
  std::string pathtofile3 = "./work/premtraction.out";
  auto file3 = std::ofstream(pathtofile3);
  for (int i = 0; i < maxmodes; ++i) {
    file3 << std::setprecision(16) << i + 1 << ";"
          << std::abs(vec_traction_r[i]) << ";"
          << std::abs(vec_traction_r_pert[i]) << ";"
          << std::abs(vec_traction_r_pertexact[i]) << std::endl;
  };

  file3.close();

  // outputting result
  std::string pathtofile4 = "./work/premcheck.out";
  auto file4 = std::ofstream(pathtofile4);
  {
    auto qx = prem_summary.q().Points();

    for (int idxe = 0; idxe < prem_summary2.fullmesh().NE(); ++idxe) {
      int idxlayer = prem_summary2.fullmesh().LayerNumber(idxe);
      for (int idxq = 0; idxq < prem_summary2.npoints(); ++idxq) {
        file4 << std::setprecision(16)
              << prem_summary2.fullmesh().NodeRadius(idxe, idxq) << ";"
              << l21_prem.u(idxlayer)(
                     prem_summary2.fullmesh().NodeRadius(idxe, idxq))
              << ";" << prem_summary2.uvec(0, idxe, idxq) << std::endl;
      }
      // };
    };
  };

  file4.close();

  // outputting result
  std::string pathtofile5 = "./work/modeintegrate1.out";
  auto file5 = std::ofstream(pathtofile5);
  {

    auto inp_model = prem;
    // Try to use integrate in gauss quad
    // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(7);
    std::vector<double> qx = q.Points();
    std::vector<double> qw = q.Weights();
    int qsize = qx.size();
    auto inp_mesh = testradmesh;

    // double totint = 0.0;
    double densitynorm = 5515.0;
    double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
    //    std::cout << radiusnorm << "\n\n";
    double pi_db = 3.14159265358979;
    double bigg_db = 6.6723 * std::pow(10.0, -11.0);
    double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
    double frequencynorm = velocitynorm / radiusnorm;

    // std::cout << "\n\n";

    for (int idxe = 0; idxe < inp_mesh.NE(); ++idxe) {
      int idxlayer = inp_mesh.LayerNumber(idxe);
      double x1 = inp_mesh.ELR(idxe);
      double x2 = inp_mesh.EUR(idxe);
      // if ((idx1 + idx2) == 0) {
      //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
      // }
      double delta = (x2 - x1) / 2.0;
      double xpl = (x1 + x2) / 2.0;
      double tmpsum = 0.0;

      // std::cout << std::setprecision(16) << x1 << " " << x2 << "\n";

      for (int idxq = 0; idxq < qsize; ++idxq) {
        // std::cout << idxe << " " << idxq << "\n";
        double xscale = inp_mesh.NodeRadius(idxe, idxq);
        // if ((x1 < 0.54) && (x2 > 0.53)) {
        //    std::cout << xscale << "\n\n";
        // }
        file5 << std::setprecision(16) << xscale << ";"
              << inp_model.Density(idxlayer)(xscale) * xscale * xscale *
                     (l21_prem.u(idxlayer)(xscale) *
                          l21_prem.u(idxlayer)(xscale) +
                      l21_prem.v(idxlayer)(xscale) *
                          l21_prem.v(idxlayer)(xscale))
              << std::endl;
      }
    };
  }
  file5.close();

  // outputting result
  std::string pathtofile6 = "./work/modeintegrate2.out";
  auto file6 = std::ofstream(pathtofile6);
  {
    auto inp_model = prem;
    auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(5);
    std::vector<double> qx = q.Points();
    std::vector<double> qw = q.Weights();
    int qsize = qx.size();
    double totint = 0.0;
    double densitynorm = 5515.0;
    double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
    //    std::cout << radiusnorm << "\n\n";
    double pi_db = 3.14159265358979;
    double bigg_db = 6.6723 * std::pow(10.0, -11.0);
    double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
    double frequencynorm = velocitynorm / radiusnorm;

    for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
      auto vec_radii = inp_model.LayerRadii(idxlayer);

      for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];
        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double tmpsum = 0.0;
        for (int idxq = 0; idxq < qsize; ++idxq) {
          double xscale = delta * qx[idxq] + xpl;

          // double tmp1 = mode1.u(idxlayer)(xscale);
          file6 << std::setprecision(16) << xscale << ";"
                << inp_model.Density(idxlayer)(xscale) *
                       (l21_prem.u(idxlayer)(xscale) *
                            l21_prem.u(idxlayer)(xscale) +
                        l21_prem.v(idxlayer)(xscale) *
                            l21_prem.v(idxlayer)(xscale)) *
                       xscale * xscale
                << std::endl;
        }
      };
    };
  }
  file6.close();

  return 0;
}