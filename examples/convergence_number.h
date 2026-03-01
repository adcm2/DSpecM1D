#include <Eigen/Core>
#include <GaussQuad/All>
#include <new_coupling/All>
#include <Interpolation/Lagrange>
#include <PlanetaryModel/All>
#include <iostream>
#include <vector>

auto
modescalarproduct(ModeCoupling::mineos_eigenfunction_continuous &mode1,
                  ModeCoupling::mineos_eigenfunction_continuous &mode2,
                  EarthModels::ModelInput<double, int> &inp_model,
                  int npoints = 5) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoints);
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

  // std::string pathtofile6 = "./work/modeintegrate2.out";
  // auto file6 = std::ofstream(pathtofile6);
  for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
    auto vec_radii = inp_model.LayerRadii(idxlayer);

    for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {
      //  if (idxlayer == inp_model.NumberOfLayers() - 1) {
      //     std::cout << vec_radii[idxpl] << " " << vec_radii[idxpl + 1]
      //               << "\n";
      //  }

      // auto fun = [&inp_model, &mode1, &mode2, &idxlayer, &vec_radii,
      //             &idxpl](double x) {
      //    double x1 = vec_radii[idxpl];
      //    double x2 = vec_radii[idxpl + 1];

      //    double delta = (x2 - x1) / 2.0;
      //    double xpl = (x1 + x2) / 2.0;
      //    double xscale = delta * x + xpl;
      //    return delta * inp_model.Density(idxlayer)(xscale) * mode1.w() *
      //           mode2.w() *
      //           (mode1.u(idxlayer)(xscale) * mode2.u(idxlayer)(xscale) +
      //            mode1.v(idxlayer)(xscale) * mode2.v(idxlayer)(xscale)) *
      //           xscale * xscale;
      // };
      double x1 = vec_radii[idxpl];
      double x2 = vec_radii[idxpl + 1];
      double delta = (x2 - x1) / 2.0;
      double xpl = (x1 + x2) / 2.0;
      double tmpsum = 0.0;

      // std::vector<double> vec_outval(101,0.0);

      // for (int idxout = 0; idxout < 1001; ++idxout) {
      //    double xscale = x1 + (x2 - x1) * idxout / 1000.0;
      //    file6 << std::setprecision(16) << xscale << ";"
      //          << inp_model.Density(idxlayer)(xscale) *
      //                 (mode1.u(idxlayer)(xscale) *
      //                      mode2.u(idxlayer)(xscale) +
      //                  mode1.v(idxlayer)(xscale) *
      //                      mode2.v(idxlayer)(xscale)) *
      //                 xscale * xscale
      //          << std::endl;
      // }
      for (int idxq = 0; idxq < qsize; ++idxq) {
        double xscale = delta * qx[idxq] + xpl;
        // double tmp1 = mode1.u(idxlayer)(xscale);

        tmpsum += qw[idxq] * inp_model.Density(idxlayer)(xscale) *
                  (mode1.u(idxlayer)(xscale) * mode2.u(idxlayer)(xscale) +
                   mode1.v(idxlayer)(xscale) * mode2.v(idxlayer)(xscale)) *
                  xscale * xscale;
      }
      // tmpsum *= delta * mode1.w() * mode2.w();
      tmpsum *= delta;
      // totint += q.Integrate(fun);
      totint += tmpsum;
    };
  };
  // file6.close();
  std::cout << frequencynorm * frequencynorm * 4.0 * pi_db * pi_db *
                   inp_model.DensityNorm() / densitynorm
            << "\n\n";
  totint *= frequencynorm * frequencynorm;
  totint *= 4.0 * pi_db * pi_db;
  totint *= inp_model.DensityNorm() / densitynorm;
  return totint;
};

auto
modescalarproduct_2(ModeCoupling::mineos_eigenfunction_continuous &mode1,
                    ModeCoupling::mineos_eigenfunction_continuous &mode2,
                    EarthModels::ModelInput<double, int> &inp_model,
                    EarthMesh::RadialMesh &inp_mesh, int npoints) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(npoints);
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

  // std::string pathtofile5 = "./work/modeintegrate1.out";
  // auto file5 = std::ofstream(pathtofile5);

  for (int idxe = 0; idxe < inp_mesh.NE(); ++idxe) {
    int idxlayer = inp_mesh.LayerNumber(idxe);
    double x1 = inp_mesh.ELR(idxe);
    double x2 = inp_mesh.EUR(idxe);
    // if ((idx1 + idx2) == 0) {
    //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
    // }
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;

    // for (int idxq = 0; idxq < qsize; ++idxq) {

    double tmpsum = 0.0;
    for (int idxq = 0; idxq < qsize; ++idxq) {
      // std::cout << idxe << " " << idxq << "\n";
      double xscale = delta * qx[idxq] + xpl;
      // double xscale = inp_mesh.NodeRadius(idxe, idxq);
      // if (idxe == 5) {
      //    std::cout << std::setprecision(16) << x1 << " " << xscale << " "
      //              << x2 << "\n";
      // }
      // file5 << std::setprecision(16) << xscale << ";"
      //       << inp_model.Density(idxlayer)(xscale) * xscale * xscale *
      //              (mode1.u(idxlayer)(xscale) * mode2.u(idxlayer)(xscale)
      //              +
      //               mode1.v(idxlayer)(xscale) * mode2.v(idxlayer)(xscale))
      //       << std::endl;
      tmpsum += qw[idxq] * inp_model.Density(idxlayer)(xscale) * xscale *
                xscale *
                (mode1.u(idxlayer)(xscale) * mode2.u(idxlayer)(xscale) +
                 mode1.v(idxlayer)(xscale) * mode2.v(idxlayer)(xscale));
    }
    tmpsum *= delta;
    totint += tmpsum;
  };
  // file5.close();
  totint *= frequencynorm * frequencynorm;
  totint *= 4.0 * pi_db * pi_db;
  totint *= inp_model.DensityNorm() / densitynorm;
  return totint;
};

auto
mode_force_vector(ModeCoupling::modecataloguecontinuous &mode_cat,
                  ModeCoupling::mineos_eigenfunction_continuous &mode_force,
                  EarthModels::ModelInput<double, int> &inp_model,
                  int maxmodes) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(7);
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();
  int qsize = qx.size();

  int maxnum = 0;
  if (maxmodes > mode_cat.NumberOfModes()) {
    maxnum = mode_cat.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  // Eigen::MatrixXd mat_gram(maxnum, maxnum);
  Eigen::VectorXd vec_force(maxnum);
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    // for (int idx2 = 0; idx2 < maxnum; ++idx2) {
    double totint = 0.0;
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
          tmpsum += qw[idxq] * inp_model.Density(idxlayer)(xscale) *
                    (mode_cat.singlemodep(idx1).u(idxlayer)(xscale) *
                         mode_force.u(idxlayer)(xscale) +
                     mode_cat.singlemodep(idx1).v(idxlayer)(xscale) *
                         mode_force.v(idxlayer)(xscale)) *
                    xscale * xscale;
          // if (((idx1 + idxlayer + idxpl) == 0) ||
          //     ((idx1 == 5) && (idxlayer == 3) && (idxpl == 0))) {
          //    std::cout
          //        << std::setprecision(16)
          //        << qw[idxq] * inp_model.Density(idxlayer)(xscale) *
          //               (mode_cat.singlemodep(idx1).u(idxlayer)(xscale)
          //               *
          //                    mode_force.u(idxlayer)(xscale) +
          //                mode_cat.singlemodep(idx1).v(idxlayer)(xscale)
          //                *
          //                    mode_force.v(idxlayer)(xscale)) *
          //               xscale * xscale
          //        << "\n";
          // }
        }
        // tmpsum *= delta * mode_cat.singlemodep(idx1).w() *
        // mode_force.w();
        tmpsum *= delta;
        // totint += q.Integrate(fun);
        totint += tmpsum;
      };
    };
    totint *= frequencynorm * frequencynorm;
    totint *= 4.0 * pi_db * pi_db;
    totint *= inp_model.DensityNorm() / densitynorm;
    vec_force(idx1) = totint;
    // }
  }
  return vec_force;
};

auto
mode_force_vector_s(ModeCoupling::modesummary &modes_all,
                    ModeCoupling::mineos_eigenfunction_continuous &mode_force,
                    EarthModels::ModelInput<double, int> &inp_model,
                    int maxmodes, bool laginc = false) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  auto q = modes_all.q();
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();
  int qsize = qx.size();

  int maxnum = 0;
  if (maxmodes > modes_all.NumberOfModes()) {
    maxnum = modes_all.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  int matsize;
  if (laginc) {
    matsize = maxnum + 2 * qsize;
  } else {
    matsize = maxnum;
  }
  // Eigen::MatrixXd mat_gram(maxnum, maxnum);
  Eigen::VectorXd vec_force(matsize);

  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    double totint = 0.0;
    for (int idxlayer = 0; idxlayer < modes_all.NumberOfLayers(); ++idxlayer) {
      auto vec_radii = inp_model.LayerRadii(idxlayer);

      for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];
        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double tmpsum = 0.0;
        for (int idxq = 0; idxq < qsize; ++idxq) {
          double xscale = delta * qx[idxq] + xpl;
          tmpsum += modes_all.multvec(idxlayer, idxpl, idxq) *
                    (modes_all.uvec(idx1, idxlayer, idxpl, idxq) *
                         mode_force.u(idxlayer)(xscale) +
                     modes_all.vvec(idx1, idxlayer, idxpl, idxq) *
                         mode_force.v(idxlayer)(xscale));
        }
        // tmpsum *= delta * modes_all.w(idx1) * modes_all.w(idx2);
        tmpsum *= delta;
        // totint += q.Integrate(fun);
        totint += tmpsum;
      };
    };
    totint *= frequencynorm * frequencynorm;
    totint *= 4.0 * pi_db * pi_db;
    totint *= inp_model.DensityNorm() / densitynorm;
    vec_force(idx1) = totint;
  }

  if (laginc) {
    int laynum = inp_model.NumberOfLayers() - 1;
    auto vec_radii = inp_model.LayerRadii(laynum);
    int lenvec = vec_radii.size();
    double deltaval = vec_radii[lenvec - 1] - vec_radii[lenvec - 2];

    // add in integral to retval
    double x1 = vec_radii[lenvec - 2];
    double x2 = vec_radii[lenvec - 1];
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;

    for (int idxl = 0; idxl < qsize; ++idxl) {
      double multval = modes_all.multvec(laynum, lenvec - 2, idxl);
      double xscale1 = delta * qx[idxl] + xpl;
      double retval1 = multval * mode_force.u(laynum)(xscale1);
      retval1 *= 2.0 * pi_db * pi_db * deltaval;
      double retval2 = multval * mode_force.v(laynum)(xscale1);
      retval2 *= 2.0 * pi_db * pi_db * deltaval;
      // if (idxl == 0) {
      //    retval1 = 0.0;
      //    retval2 = 0.0;
      // }
      vec_force(maxnum + idxl) = retval1;
      vec_force(maxnum + qsize + idxl) = retval2;
    }
  }

  return vec_force;
};

auto
mode_force_vector_s2(ModeCoupling::modesummary2 &modes_all,
                     ModeCoupling::mineos_eigenfunction_continuous &mode_force,
                     EarthModels::ModelInput<double, int> &inp_model,
                     int maxmodes, bool laginc = false) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  // auto q = modes_all.q();
  // std::vector<double> qx = q.Points();
  // std::vector<double> qw = q.Weights();
  int qsize = modes_all.npoints();

  int maxnum = 0;
  if (maxmodes > modes_all.NumberOfModes()) {
    maxnum = modes_all.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  int matsize;
  if (laginc) {
    matsize = maxnum + 2 * qsize;
  } else {
    matsize = maxnum;
  }
  // Eigen::MatrixXd mat_gram(maxnum, maxnum);
  Eigen::VectorXd vec_force(matsize);

  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    double totint = 0.0;

    // double totint = 0.0;

    for (int idxe = 0; idxe < modes_all.fullmesh().NE(); ++idxe) {
      int idxlayer = modes_all.fullmesh().LayerNumber(idxe);
      double x1 = modes_all.fullmesh().ELR(idxe);
      double x2 = modes_all.fullmesh().EUR(idxe);
      // if ((idx1 + idx2) == 0) {
      //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
      // }
      double delta = (x2 - x1) / 2.0;
      double xpl = (x1 + x2) / 2.0;
      double tmpsum = 0.0;
      for (int idxq = 0; idxq < qsize; ++idxq) {
        // std::cout << idxe << " " << idxq << "\n";
        double xscale = modes_all.fullmesh().NodeRadius(idxe, idxq);
        tmpsum +=
            modes_all.multvec(idxe, idxq) *
            (modes_all.uvec(idx1, idxe, idxq) * mode_force.u(idxlayer)(xscale) +
             modes_all.vvec(idx1, idxe, idxq) * mode_force.v(idxlayer)(xscale));
      }
      tmpsum *= delta;
      totint += tmpsum;
      // if (idxe < modes_all.fullmesh().NE() - 1) {
      //    if (((idx1 + idx2) == 0) && (modes_all.fullmesh().LayerNumber(
      //                                     idxe + 1) == (idxlayer + 1)))
      //                                     {
      //       std::cout << std::setprecision(16) << totint << "\n";
      //    }
      // }
      // };
    };
    totint *= frequencynorm * frequencynorm;
    totint *= 4.0 * pi_db * pi_db;
    totint *= inp_model.DensityNorm() / densitynorm;
    vec_force(idx1) = totint;
  }

  if (laginc) {
    int laynum = inp_model.NumberOfLayers() - 1;
    auto vec_radii = inp_model.LayerRadii(laynum);
    int lenvec = vec_radii.size();
    double deltaval = vec_radii[lenvec - 1] - vec_radii[lenvec - 2];

    // add in integral to retval
    double x1 = vec_radii[lenvec - 2];
    double x2 = vec_radii[lenvec - 1];
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;

    for (int idxl = 0; idxl < qsize; ++idxl) {
      // double multval = modes_all.multvec(laynum, lenvec - 2, idxl);
      // double xscale1 = delta * qx[idxl] + xpl;
      // double retval1 = multval * mode_force.u(laynum)(xscale1);
      // retval1 *= 2.0 * pi_db * pi_db * deltaval;
      // double retval2 = multval * mode_force.v(laynum)(xscale1);
      // retval2 *= 2.0 * pi_db * pi_db * deltaval;
      // // if (idxl == 0) {
      // //    retval1 = 0.0;
      // //    retval2 = 0.0;
      // // }
      // vec_force(maxnum + idxl) = retval1;
      // vec_force(maxnum + qsize + idxl) = retval2;
    }
  }

  if (laginc) {
    int elemnum = modes_all.fullmesh().NE() - 1;
    // int laynum = inp_model.NumberOfLayers() - 1;
    // auto vec_radii = inp_model.LayerRadii(laynum);
    // int lenvec = vec_radii.size();
    double deltaval = modes_all.fullmesh().EW(elemnum);
    int laynum = modes_all.fullmesh().LayerNumber(elemnum);
    // add in integral to retval
    // double x1 = modes_all.fullmesh().ELR(elemnum);
    // double x2 = modes_all.fullmesh().EUR(elemnum);
    // double delta = (x2 - x1) / 2.0;
    // double xpl = (x1 + x2) / 2.0;

    // mode against lagrange polynomial
    for (int idxl = 0; idxl < qsize; ++idxl) {
      // double xscale1 = delta * qx[idxl] + xpl;
      // double xscale2 = delta * qx[idxl] + xpl;
      double multval = modes_all.multvec(elemnum, idxl);
      // double wdr = qw[idxl] * inp_model.Density(laynum)(xscale1)
      // for (int idxm = 0; idxm < maxnum; ++idxm) {
      double retval1 =
          multval *
          mode_force.u(laynum)(modes_all.fullmesh().NodeRadius(elemnum, idxl));
      retval1 *= 2.0 * pi_db * pi_db * deltaval;
      vec_force(maxnum + idxl) = retval1;

      double retval2 =
          multval *
          mode_force.v(laynum)(modes_all.fullmesh().NodeRadius(elemnum, idxl));
      retval2 *= 2.0 * pi_db * pi_db * deltaval;
      vec_force(maxnum + qsize + idxl) = retval2;
      // }
    }
  }

  return vec_force;
};

auto
mode_force_vector_s3(ModeCoupling::modesummary2 &modes_all,
                     ModeCoupling::mineos_eigenfunction_continuous &mode_force,
                     EarthModels::ModelInput<double, int> &inp_model,
                     int maxmodes, Eigen::VectorXd &vec_coeff) {
  int qsize = modes_all.npoints();

  int maxnum = 0;
  if (maxmodes > modes_all.NumberOfModes()) {
    maxnum = modes_all.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  int matsize;
  matsize = maxnum;
  Eigen::VectorXd vec_force(matsize);

  // assumption that there are only 4 augmented functions and that they are the
  // sin(x) and sin(2x) used in the derivations
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    double totint = 0.0;

    for (int idxe = 0; idxe < modes_all.fullmesh().NE(); ++idxe) {
      int idxlayer = modes_all.fullmesh().LayerNumber(idxe);

      double x1 = modes_all.fullmesh().ELR(idxe);
      double x2 = modes_all.fullmesh().EUR(idxe);
      double delta = (x2 - x1) / 2.0;
      double xpl = (x1 + x2) / 2.0;
      double tmpsum = 0.0;
      for (int idxq = 0; idxq < qsize; ++idxq) {
        // std::cout << idxe << " " << idxq << "\n";
        double xscale = modes_all.fullmesh().NodeRadius(idxe, idxq);
        double tmpu = 0.0, tmpv = 0.0;
        if (idxlayer == inp_model.NumberOfLayers() - 1) {
          double xlow = inp_model.LowerRadius(idxlayer);
          double xup = inp_model.UpperRadius(idxlayer);
          tmpu +=
              vec_coeff(0) * std::sin(pi_db * (xscale - xlow) / (xup - xlow));
          tmpu += vec_coeff(1) *
                  std::sin(2.0 * pi_db * (xscale - xlow) / (xup - xlow));
          tmpv +=
              vec_coeff(2) * std::sin(pi_db * (xscale - xlow) / (xup - xlow));
          tmpv += vec_coeff(3) *
                  std::sin(2.0 * pi_db * (xscale - xlow) / (xup - xlow));
        }

        tmpsum += modes_all.multvec(idxe, idxq) *
                  (modes_all.uvec(idx1, idxe, idxq) *
                       (mode_force.u(idxlayer)(xscale) - tmpu) +
                   modes_all.vvec(idx1, idxe, idxq) *
                       (mode_force.v(idxlayer)(xscale) - tmpv));
      }
      tmpsum *= delta;
      totint += tmpsum;
    };
    totint *= frequencynorm * frequencynorm;
    totint *= 4.0 * pi_db * pi_db;
    totint *= inp_model.DensityNorm() / densitynorm;
    vec_force(idx1) = totint;
  }

  return vec_force;
};

auto
mode_force_lagrange_vector(
    ModeCoupling::mineos_eigenfunction_continuous &mode_force,
    EarthModels::ModelInput<double, int> &inp_model, int maxmodes,
    std::vector<int> lagint) {

  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(7);
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();
  int qsize = qx.size();

  // int maxnum = 0;
  // if (maxmodes > mode_cat.NumberOfModes()) {
  //    maxnum = mode_cat.NumberOfModes();
  // } else {
  //    maxnum = maxmodes;
  // }

  // double densitynorm = 5515.0;
  // double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  // //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  // double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  // double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db *
  // densitynorm); double frequencynorm = velocitynorm / radiusnorm;

  // // Eigen::MatrixXd mat_gram(maxnum, maxnum);
  // Eigen::VectorXd vec_force(maxnum);

  // for (int idx1 = 0; idx1 < maxnum; ++idx1) {
  //    // for (int idx2 = 0; idx2 < maxnum; ++idx2) {
  //    double totint = 0.0;
  int laynum = inp_model.NumberOfLayers() - 1;
  auto vec_radii = inp_model.LayerRadii(laynum);
  int lenvec = vec_radii.size();
  double deltaval = vec_radii[lenvec - 1] - vec_radii[lenvec - 2];

  // add in integral to retval
  double x1 = vec_radii[lenvec - 2];
  double x2 = vec_radii[lenvec - 1];
  double delta = (x2 - x1) / 2.0;
  double xpl = (x1 + x2) / 2.0;
  double xscale1 = delta * qx[lagint[0]] + xpl;
  double xscale2 = delta * qx[lagint[1]] + xpl;
  double retval = 0.0;
  retval += qw[lagint[0]] * inp_model.Density(laynum)(xscale1) *
            mode_force.u(laynum)(xscale1) * xscale1 * xscale1;
  retval += qw[lagint[1]] * inp_model.Density(laynum)(xscale2) *
            mode_force.v(laynum)(xscale2) * xscale2 * xscale2;
  retval *= 2.0 * pi_db * pi_db * deltaval;

  return retval;
};

auto
mode_gram_matrix(ModeCoupling::modecataloguecontinuous &mode_cat,
                 EarthModels::ModelInput<double, int> &inp_model, int maxmodes,
                 bool laginc = false) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();
  int qsize = qx.size();

  int maxnum = 0;
  if (maxmodes > mode_cat.NumberOfModes()) {
    maxnum = mode_cat.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  //    if (laginc) {
  //       std::cout << "Included\n\n";
  //    } else {
  //       std::cout << "Not included\n\n";
  //    }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  Eigen::MatrixXd mat_gram(maxnum, maxnum);

  // trying slightly quicker method
  // std::vector<std::vector<std::vector<double>>> vec_wdrr(
  //     inp_model.NumberOfLayers(),
  //     std::vector<std::vector<double>>(inp_model.LayerRadii(0).size(),
  //                                      std::vector<double>(qsize, 0.0)));
  std::vector<std::vector<std::vector<double>>> vec_wdrr;
  for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
    auto vec_radii = inp_model.LayerRadii(idxlayer);
    std::vector<std::vector<double>> vec_tmp1;
    // std::cout << idxlayer << "\n";
    for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

      double x1 = vec_radii[idxpl];
      double x2 = vec_radii[idxpl + 1];
      double delta = (x2 - x1) / 2.0;
      double xpl = (x1 + x2) / 2.0;
      double tmpsum = 0.0;
      std::vector<double> vec_tmp2(qsize, 0.0);
      for (int idxq = 0; idxq < qsize; ++idxq) {
        double xscale = delta * qx[idxq] + xpl;
        // double tmp1 = mode1.u(idxlayer)(xscale);
        // tmpsum += qw[idxq] * inp_model.Density(idxlayer)(xscale) *
        //           (mode_cat.singlemodep(idx1).u(idxlayer)(xscale) *
        //                mode_cat.singlemodep(idx2).u(idxlayer)(xscale) +
        //            mode_cat.singlemodep(idx1).v(idxlayer)(xscale) *
        //                mode_cat.singlemodep(idx2).v(idxlayer)(xscale)) *
        //           xscale * xscale;
        vec_tmp2[idxq] =
            qw[idxq] * inp_model.Density(idxlayer)(xscale) * xscale * xscale;
      }
      vec_tmp1.push_back(vec_tmp2);
      // tmpsum *= delta * mode_cat.singlemodep(idx1).w() *
      //           mode_cat.singlemodep(idx2).w();
      // // totint += q.Integrate(fun);
      // totint += tmpsum;
    };
    vec_wdrr.push_back(vec_tmp1);
  };
  std::vector<std::vector<std::vector<std::vector<double>>>> vec_uval, vec_vval;
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    std::vector<std::vector<std::vector<double>>> vec_tmp1, vec_tmp11;
    for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers(); ++idxlayer) {
      auto vec_radii = inp_model.LayerRadii(idxlayer);
      std::vector<std::vector<double>> vec_tmp2, vec_tmp21;
      for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];
        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        // double tmpsum = 0.0;
        std::vector<double> vec_tmp3(qsize, 0.0), vec_tmp31(qsize, 0.0);
        for (int idxq = 0; idxq < qsize; ++idxq) {
          double xscale = delta * qx[idxq] + xpl;
          // double tmp1 = mode1.u(idxlayer)(xscale);
          // tmpsum +=
          //     qw[idxq] * inp_model.Density(idxlayer)(xscale) *
          //     (mode_cat.singlemodep(idx1).u(idxlayer)(xscale) *
          //          mode_cat.singlemodep(idx2).u(idxlayer)(xscale) +
          //      mode_cat.singlemodep(idx1).v(idxlayer)(xscale) *
          //          mode_cat.singlemodep(idx2).v(idxlayer)(xscale)) *
          //     xscale * xscale;
          vec_tmp3[idxq] = mode_cat.singlemodep(idx1).u(idxlayer)(xscale);
          vec_tmp31[idxq] = mode_cat.singlemodep(idx1).v(idxlayer)(xscale);
        }
        vec_tmp2.push_back(vec_tmp3);
        vec_tmp21.push_back(vec_tmp31);
      };
      vec_tmp1.push_back(vec_tmp2);
      vec_tmp11.push_back(vec_tmp21);
    };
    vec_uval.push_back(vec_tmp1);
    vec_vval.push_back(vec_tmp11);
  }
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    for (int idx2 = 0; idx2 < maxnum; ++idx2) {
      double totint = 0.0;
      for (int idxlayer = 0; idxlayer < inp_model.NumberOfLayers();
           ++idxlayer) {
        auto vec_radii = inp_model.LayerRadii(idxlayer);

        for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

          double x1 = vec_radii[idxpl];
          double x2 = vec_radii[idxpl + 1];
          double delta = (x2 - x1) / 2.0;
          double xpl = (x1 + x2) / 2.0;
          double tmpsum = 0.0;
          for (int idxq = 0; idxq < qsize; ++idxq) {
            double xscale = delta * qx[idxq] + xpl;
            tmpsum += vec_wdrr[idxlayer][idxpl][idxq] *
                      (vec_uval[idx1][idxlayer][idxpl][idxq] *
                           vec_uval[idx2][idxlayer][idxpl][idxq] +
                       vec_vval[idx1][idxlayer][idxpl][idxq] *
                           vec_vval[idx2][idxlayer][idxpl][idxq]);
          }
          // tmpsum *= delta * mode_cat.singlemodep(idx1).w() *
          //           mode_cat.singlemodep(idx2).w();
          tmpsum *= delta;
          // totint += q.Integrate(fun);
          totint += tmpsum;
        };
      };
      totint *= frequencynorm * frequencynorm;
      totint *= 4.0 * pi_db * pi_db;
      totint *= inp_model.DensityNorm() / densitynorm;
      mat_gram(idx1, idx2) = totint;
    }
  }
  return mat_gram;
};

auto
mode_gram_matrix_s(ModeCoupling::modesummary &modes_all,
                   EarthModels::ModelInput<double, int> &inp_model,
                   int maxmodes, bool laginc = false) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  auto q = modes_all.q();
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();
  int qsize = modes_all.npoints();

  int maxnum = 0;
  if (maxmodes > modes_all.NumberOfModes()) {
    maxnum = modes_all.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  //    if (laginc) {
  //       std::cout << "Included\n\n";
  //    } else {
  //       std::cout << "Not included\n\n";
  //    }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  int matsize;
  if (laginc) {
    matsize = maxnum + 2 * qsize;
  } else {
    matsize = maxnum;
  }

  Eigen::MatrixXd mat_gram(matsize, matsize);
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    for (int idx2 = 0; idx2 < maxnum; ++idx2) {
      double totint = 0.0;
      // int numlaytot = 0;
      for (int idxlayer = 0; idxlayer < modes_all.NumberOfLayers();
           ++idxlayer) {
        auto vec_radii = inp_model.LayerRadii(idxlayer);

        for (int idxpl = 0; idxpl < vec_radii.size() - 1; ++idxpl) {

          double x1 = vec_radii[idxpl];
          double x2 = vec_radii[idxpl + 1];
          double delta = (x2 - x1) / 2.0;
          double xpl = (x1 + x2) / 2.0;
          double tmpsum = 0.0;
          // if ((idx1 + idx2) == 0) {
          //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
          // }
          for (int idxq = 0; idxq < qsize; ++idxq) {
            tmpsum += modes_all.multvec(idxlayer, idxpl, idxq) *
                      (modes_all.uvec(idx1, idxlayer, idxpl, idxq) *
                           modes_all.uvec(idx2, idxlayer, idxpl, idxq) +
                       modes_all.vvec(idx1, idxlayer, idxpl, idxq) *
                           modes_all.vvec(idx2, idxlayer, idxpl, idxq));
          }
          // tmpsum *= delta * modes_all.w(idx1) * modes_all.w(idx2);
          tmpsum *= delta;
          // totint += q.Integrate(fun);
          totint += tmpsum;
          // ++numlaytot;
        };
        // if ((idx1 + idx2) == 0) {
        //    std::cout << std::setprecision(16) << totint << "\n";
        // }
      };
      totint *= frequencynorm * frequencynorm;
      totint *= 4.0 * pi_db * pi_db;
      totint *= inp_model.DensityNorm() / densitynorm;
      mat_gram(idx1, idx2) = totint;

      // if ((idx1 + idx2) == 0) {
      //    std::cout << "Number of elements: " << numlaytot << "\n";
      // }
    }
  }

  // if we are including the lagrange polynomials
  if (laginc) {
    int laynum = inp_model.NumberOfLayers() - 1;
    auto vec_radii = inp_model.LayerRadii(laynum);
    int lenvec = vec_radii.size();
    double deltaval = vec_radii[lenvec - 1] - vec_radii[lenvec - 2];

    // add in integral to retval
    double x1 = vec_radii[lenvec - 2];
    double x2 = vec_radii[lenvec - 1];
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;
    // mode against lagrange polynomial
    for (int idxl = 0; idxl < qsize; ++idxl) {
      double xscale1 = delta * qx[idxl] + xpl;
      double xscale2 = delta * qx[idxl] + xpl;
      double multval = modes_all.multvec(laynum, lenvec - 2, idxl);
      // double wdr = qw[idxl] * inp_model.Density(laynum)(xscale1)
      for (int idxm = 0; idxm < maxnum; ++idxm) {
        double retval1 =
            multval * modes_all.uvec(idxm, laynum, lenvec - 2, idxl);
        retval1 *= 2.0 * pi_db * pi_db * deltaval;
        mat_gram(maxnum + idxl, idxm) = retval1;
        mat_gram(idxm, maxnum + idxl) = retval1;

        double retval2 =
            multval * modes_all.vvec(idxm, laynum, lenvec - 2, idxl);
        retval2 *= 2.0 * pi_db * pi_db * deltaval;
        mat_gram(maxnum + qsize + idxl, idxm) = retval2;
        mat_gram(idxm, maxnum + idxl + qsize) = retval2;
      }
    }

    // lagrange polynomial against each other
    for (int idxl = 0; idxl < qsize; ++idxl) {
      double retval = modes_all.multvec(laynum, lenvec - 2, idxl);
      retval *= 2.0 * pi_db * pi_db * deltaval;
      //  std::cout << idxl << " " << retval << "\n";
      mat_gram(maxnum + idxl, maxnum + idxl) = retval;
      mat_gram(maxnum + qsize + idxl, maxnum + qsize + idxl) = retval;
    }
  }
  return mat_gram;
};

auto
mode_gram_matrix_s2(ModeCoupling::modesummary2 &modes_all,
                    EarthModels::ModelInput<double, int> &inp_model,
                    int maxmodes, bool laginc = false) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);
  // auto q = modes_all.q();
  // std::vector<double> qx = q.Points();
  // std::vector<double> qw = q.Weights();
  int qsize = modes_all.npoints();

  int maxnum = 0;
  if (maxmodes > modes_all.NumberOfModes()) {
    maxnum = modes_all.NumberOfModes();
  } else {
    maxnum = maxmodes;
  }

  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  int matsize;
  if (laginc) {
    matsize = maxnum + 2 * qsize;
  } else {
    matsize = maxnum;
  }
  Eigen::MatrixXd mat_gram(matsize, matsize);
  for (int idx1 = 0; idx1 < maxnum; ++idx1) {
    for (int idx2 = 0; idx2 < maxnum; ++idx2) {
      // std::cout << "idx1/2: " << idx1 << " " << idx2 << "\n";
      double totint = 0.0;

      for (int idxe = 0; idxe < modes_all.fullmesh().NE(); ++idxe) {
        int idxlayer = modes_all.fullmesh().LayerNumber(idxe);
        double x1 = modes_all.fullmesh().ELR(idxe);
        double x2 = modes_all.fullmesh().EUR(idxe);
        // if ((idx1 + idx2) == 0) {
        //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
        // }
        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double tmpsum = 0.0;
        for (int idxq = 0; idxq < qsize; ++idxq) {
          // std::cout << idxe << " " << idxq << "\n";
          tmpsum += modes_all.multvec(idxe, idxq) *
                    (modes_all.uvec(idx1, idxe, idxq) *
                         modes_all.uvec(idx2, idxe, idxq) +
                     modes_all.vvec(idx1, idxe, idxq) *
                         modes_all.vvec(idx2, idxe, idxq));
        }
        tmpsum *= delta;
        totint += tmpsum;
      };
      totint *= frequencynorm * frequencynorm;
      totint *= 4.0 * pi_db * pi_db;
      totint *= inp_model.DensityNorm() / densitynorm;
      mat_gram(idx1, idx2) = totint;
    }
  }

  // if we are including the lagrange polynomials
  if (laginc) {
    int elemnum = modes_all.fullmesh().NE() - 1;
    double deltaval = modes_all.fullmesh().EW(elemnum);

    // add in integral to retval
    double x1 = modes_all.fullmesh().ELR(elemnum);
    double x2 = modes_all.fullmesh().EUR(elemnum);
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;

    // mode against lagrange polynomial
    for (int idxl = 0; idxl < qsize; ++idxl) {
      double multval = modes_all.multvec(elemnum, idxl);
      for (int idxm = 0; idxm < maxnum; ++idxm) {
        double retval1 = multval * modes_all.uvec(idxm, elemnum, idxl);
        retval1 *= 2.0 * pi_db * pi_db * deltaval;
        mat_gram(maxnum + idxl, idxm) = retval1;
        mat_gram(idxm, maxnum + idxl) = retval1;

        double retval2 = multval * modes_all.vvec(idxm, elemnum, idxl);
        retval2 *= 2.0 * pi_db * pi_db * deltaval;
        mat_gram(maxnum + qsize + idxl, idxm) = retval2;
        mat_gram(idxm, maxnum + idxl + qsize) = retval2;
      }
    }

    // lagrange polynomial against each other
    for (int idxl = 0; idxl < qsize; ++idxl) {
      double retval = modes_all.multvec(elemnum, idxl);
      retval *= 2.0 * pi_db * pi_db * deltaval;
      mat_gram(maxnum + idxl, maxnum + idxl) = retval;
      mat_gram(maxnum + qsize + idxl, maxnum + qsize + idxl) = retval;
    }
  }
  return mat_gram;
};

auto
l2norm(ModeCoupling::mineos_eigenfunction_continuous &mode1,
       EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();

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
      auto fun = [&inp_model, &mode1, &idxlayer, &vec_radii, &idxpl](double x) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];

        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double xscale = delta * x + xpl;
        return delta * mode1.u(idxlayer)(xscale) * mode1.u(idxlayer)(xscale);
      };
      totint += q.Integrate(fun);
    };
  };
  return std::sqrt(totint);
};

auto
l2norm_v(ModeCoupling::mineos_eigenfunction_continuous &mode1,
         EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();

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
      auto fun = [&inp_model, &mode1, &idxlayer, &vec_radii, &idxpl](double x) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];

        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double xscale = delta * x + xpl;
        return delta * mode1.v(idxlayer)(xscale) * mode1.v(idxlayer)(xscale);
      };
      totint += q.Integrate(fun);
    };
  };
  return std::sqrt(totint);
};

auto
modesum(Eigen::VectorXd &mode_coefficients,
        ModeCoupling::modecataloguecontinuous &mode_cat,
        EarthModels::ModelInput<double, int> &inp_model, int laynum,
        double radius) {
  double sumval = 0.0;
  for (int idx = 0; idx < mode_coefficients.rows(); ++idx) {
    // double tmp = mode_cat.singlemode(idx).u(idxlayer)(xscale);
    sumval +=
        mode_coefficients(idx) * mode_cat.singlemodep(idx).u(laynum)(radius);
  }

  return sumval;
};

auto
modesum_s(Eigen::VectorXd &mode_coefficients,
          ModeCoupling::modesummary &mode_sum,
          EarthModels::ModelInput<double, int> &inp_model, int laynum,
          int idxpl, int idxq, bool laginc = false) {
  double sumval = 0.0;
  // auto q = mode_sum.q();
  int qsize = mode_sum.npoints();
  int nummodes;
  if (laginc) {
    nummodes = mode_coefficients.rows() - 2 * qsize;
  } else {
    nummodes = mode_coefficients.rows();
  }

  for (int idx = 0; idx < nummodes; ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.uvec(idx, laynum, idxpl, idxq);
  }
  if (laginc) {
    if (laynum == (inp_model.NumberOfLayers() - 1)) {
      if (idxpl == (inp_model.LayerRadiiNumber(laynum) - 2)) {
        // if (idxq != 0) {
        sumval += mode_coefficients(nummodes + idxq);
        // }
      }
    }
  }

  return sumval;
};

auto
modesum_s2(Eigen::VectorXd &mode_coefficients,
           ModeCoupling::modesummary2 &mode_sum,
           EarthModels::ModelInput<double, int> &inp_model, int elemnum,
           int idxq, bool laginc = false) {
  double sumval = 0.0;
  // auto q = mode_sum.q();
  int qsize = mode_sum.npoints();
  int nummodes;
  if (laginc) {
    nummodes = mode_coefficients.rows() - 2 * qsize;
  } else {
    nummodes = mode_coefficients.rows();
  }

  for (int idx = 0; idx < nummodes; ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.uvec(idx, elemnum, idxq);
  }
  if (laginc) {
    if (elemnum == (mode_sum.fullmesh().NE() - 1)) {
      // if (idxq != 0) {
      sumval += mode_coefficients(nummodes + idxq);
    }
  }

  return sumval;
};

auto
modesum_sup(Eigen::VectorXd &mode_coefficients,
            ModeCoupling::modesummary &mode_sum,
            EarthModels::ModelInput<double, int> &inp_model, int laynum,
            int idxpl, int idxq) {
  double sumval = 0.0;
  for (int idx = 0; idx < mode_coefficients.rows(); ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.upvec(idx, laynum, idxpl, idxq);
  }

  return sumval;
};

auto
modesum_sv(Eigen::VectorXd &mode_coefficients,
           ModeCoupling::modesummary &mode_sum,
           EarthModels::ModelInput<double, int> &inp_model, int laynum,
           int idxpl, int idxq, bool laginc = false) {
  double sumval = 0.0;
  int qsize = mode_sum.npoints();
  int nummodes;
  if (laginc) {
    nummodes = mode_coefficients.rows() - 2 * qsize;
  } else {
    nummodes = mode_coefficients.rows();
  }

  for (int idx = 0; idx < nummodes; ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.vvec(idx, laynum, idxpl, idxq);
  }
  if (laginc) {
    if (laynum == (inp_model.NumberOfLayers() - 1)) {
      if (idxpl == (inp_model.LayerRadiiNumber(laynum) - 2)) {
        // if (idxq != 0) {
        sumval += mode_coefficients(nummodes + qsize + idxq);
        // }
      }
    }
  }
  return sumval;
};

auto
modesum_sv2(Eigen::VectorXd &mode_coefficients,
            ModeCoupling::modesummary2 &mode_sum,
            EarthModels::ModelInput<double, int> &inp_model, int elemnum,
            int idxq, bool laginc = false) {
  double sumval = 0.0;
  // auto q = mode_sum.q();
  int qsize = mode_sum.npoints();
  int nummodes;
  if (laginc) {
    nummodes = mode_coefficients.rows() - 2 * qsize;
  } else {
    nummodes = mode_coefficients.rows();
  }

  for (int idx = 0; idx < nummodes; ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.vvec(idx, elemnum, idxq);
  }
  if (laginc) {
    if (elemnum == (mode_sum.fullmesh().NE() - 1)) {
      // if (idxq != 0) {
      sumval += mode_coefficients(nummodes + qsize + idxq);
    }
  }

  return sumval;
};

auto
modesum_sup2(Eigen::VectorXd &mode_coefficients,
             ModeCoupling::modesummary2 &mode_sum,
             EarthModels::ModelInput<double, int> &inp_model, int elemnum,
             int idxq, bool laginc = false) {
  double sumval = 0.0;
  // auto q = mode_sum.q();
  int qsize = mode_sum.npoints();
  int nummodes;
  if (laginc) {
    nummodes = mode_coefficients.rows() - 2 * qsize;
  } else {
    nummodes = mode_coefficients.rows();
  }

  for (int idx = 0; idx < nummodes; ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.upvec(idx, elemnum, idxq);
  }
  if (laginc) {
    if (elemnum == (mode_sum.fullmesh().NE() - 1)) {
      // if (idxq != 0) {
      auto pleg = Interpolation::LagrangePolynomial(
          mode_sum.q().Points().begin(), mode_sum.q().Points().end());
      auto xscale = mode_sum.q().X(idxq);
      for (int i = 0; i < qsize; ++i) {
        sumval += mode_coefficients(nummodes + idxq) *
                  pleg.Derivative(i, xscale) * 2.0 /
                  mode_sum.fullmesh().EW(elemnum);
      }
    }
  }

  return sumval;
};

auto
modesum_svp2(Eigen::VectorXd &mode_coefficients,
             ModeCoupling::modesummary2 &mode_sum,
             EarthModels::ModelInput<double, int> &inp_model, int elemnum,
             int idxq, bool laginc = false) {
  double sumval = 0.0;
  // auto q = mode_sum.q();
  int qsize = mode_sum.npoints();
  int nummodes;
  if (laginc) {
    nummodes = mode_coefficients.rows() - 2 * qsize;
  } else {
    nummodes = mode_coefficients.rows();
  }

  for (int idx = 0; idx < nummodes; ++idx) {
    sumval += mode_coefficients(idx) * mode_sum.vpvec(idx, elemnum, idxq);
  }
  if (laginc) {
    if (elemnum == (mode_sum.fullmesh().NE() - 1)) {
      // if (idxq != 0) {
      auto pleg = Interpolation::LagrangePolynomial(
          mode_sum.q().Points().begin(), mode_sum.q().Points().end());
      auto xscale = mode_sum.q().X(idxq);
      for (int i = 0; i < qsize; ++i) {
        sumval += mode_coefficients(nummodes + qsize + idxq) *
                  pleg.Derivative(i, xscale) * 2.0 /
                  mode_sum.fullmesh().EW(elemnum);
      }
    }
  }

  return sumval;
};

auto
l2normbasissum(Eigen::VectorXd &mode_coefficients,
               ModeCoupling::modecataloguecontinuous &mode_cat,
               EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);

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
      //  if (idxlayer == inp_model.NumberOfLayers() - 1) {
      //     std::cout << vec_radii[idxpl] << " " << vec_radii[idxpl + 1]
      //               << "\n";
      //  }

      auto fun = [&inp_model, &mode_coefficients, &mode_cat, &idxlayer,
                  &vec_radii, &idxpl](double x) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];

        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double xscale = delta * x + xpl;
        double sumval = 0.0;
        double tmpval =
            modesum(mode_coefficients, mode_cat, inp_model, idxlayer, xscale);
        return delta * tmpval * tmpval;
      };
      totint += q.Integrate(fun);
    };
  };
  // totint *= frequencynorm * frequencynorm;
  // totint *= 4.0 * pi_db * pi_db;
  // totint *= inp_model.DensityNorm() / densitynorm;
  return std::sqrt(totint);
};

auto
l2normdifference(Eigen::VectorXd &mode_coefficients,
                 ModeCoupling::modecataloguecontinuous &mode_cat,
                 ModeCoupling::mineos_eigenfunction_continuous &modecomp,
                 EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(4);

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
      //  if (idxlayer == inp_model.NumberOfLayers() - 1) {
      //     std::cout << vec_radii[idxpl] << " " << vec_radii[idxpl + 1]
      //               << "\n";
      //  }

      auto fun = [&inp_model, &mode_coefficients, &mode_cat, &modecomp,
                  &idxlayer, &vec_radii, &idxpl](double x) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];

        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double xscale = delta * x + xpl;
        double sumval = 0.0;
        double tmpval =
            modesum(mode_coefficients, mode_cat, inp_model, idxlayer, xscale);
        tmpval -= modecomp.u(idxlayer)(xscale);
        return delta * tmpval * tmpval;
      };
      totint += q.Integrate(fun);
    };
  };
  // totint *= frequencynorm * frequencynorm;
  // totint *= 4.0 * pi_db * pi_db;
  // totint *= inp_model.DensityNorm() / densitynorm;
  return std::sqrt(totint);
};

auto
l2normdifference_s(Eigen::VectorXd &mode_coefficients,
                   ModeCoupling::modesummary &mode_all,
                   ModeCoupling::mineos_eigenfunction_continuous &modecomp,
                   EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = mode_all.q().Points();
  std::vector<double> qw = mode_all.q().Weights();
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
      for (int idxq = 0; idxq < qx.size(); ++idxq) {
        double tmpval = modesum_s(mode_coefficients, mode_all, inp_model,
                                  idxlayer, idxpl, idxq);
        double xscale = delta * qx[idxq] + xpl;
        tmpval -= modecomp.u(idxlayer)(xscale);
        tmpsum += delta * tmpval * tmpval * qw[idxq];
      }
      totint += tmpsum;
    };
  };
  return std::sqrt(totint);
};

auto
l2normdifference_s_full(Eigen::VectorXd &mode_coefficients,
                        ModeCoupling::modesummary &mode_all,
                        ModeCoupling::mineos_eigenfunction_continuous &modecomp,
                        EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = mode_all.q().Points();
  std::vector<double> qw = mode_all.q().Weights();
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
      for (int idxq = 0; idxq < qx.size(); ++idxq) {
        double tmpval = modesum_s(mode_coefficients, mode_all, inp_model,
                                  idxlayer, idxpl, idxq);
        double tmpval2 = modesum_sv(mode_coefficients, mode_all, inp_model,
                                    idxlayer, idxpl, idxq);
        double xscale = delta * qx[idxq] + xpl;
        tmpval -= modecomp.u(idxlayer)(xscale);
        tmpval2 -= modecomp.v(idxlayer)(xscale);
        tmpsum += delta * qw[idxq] * xscale * xscale *
                  inp_model.Density(idxlayer)(xscale) *
                  (tmpval * tmpval + tmpval2 * tmpval2);
      }
      totint += tmpsum;
    };
  };
  return std::sqrt(totint);
};

auto
l2normdifference_s_full2(
    Eigen::VectorXd &mode_coefficients, ModeCoupling::modesummary2 &modes_all,
    ModeCoupling::mineos_eigenfunction_continuous &modecomp,
    EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = modes_all.q().Points();
  std::vector<double> qw = modes_all.q().Weights();
  double totint = 0.0;
  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  for (int idxe = 0; idxe < modes_all.fullmesh().NE(); ++idxe) {
    int idxlayer = modes_all.fullmesh().LayerNumber(idxe);
    double x1 = modes_all.fullmesh().ELR(idxe);
    double x2 = modes_all.fullmesh().EUR(idxe);
    // if ((idx1 + idx2) == 0) {
    //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
    // }
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;
    double tmpsum = 0.0;
    for (int idxq = 0; idxq < qx.size(); ++idxq) {

      double tmpval =
          modesum_s2(mode_coefficients, modes_all, inp_model, idxe, idxq);
      double tmpval2 =
          modesum_sv2(mode_coefficients, modes_all, inp_model, idxe, idxq);
      double xscale = delta * qx[idxq] + xpl;
      tmpval -= modecomp.u(idxlayer)(xscale);
      tmpval2 -= modecomp.v(idxlayer)(xscale);
      tmpsum += delta * qw[idxq] * xscale * xscale *
                inp_model.Density(idxlayer)(xscale) *
                (tmpval * tmpval + tmpval2 * tmpval2);
    }
    totint += tmpsum;
  };
  return std::sqrt(totint);
};

auto
l2normdifference_sv(Eigen::VectorXd &mode_coefficients,
                    ModeCoupling::modesummary &mode_all,
                    ModeCoupling::mineos_eigenfunction_continuous &modecomp,
                    EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = mode_all.q().Points();
  std::vector<double> qw = mode_all.q().Weights();
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
      for (int idxq = 0; idxq < qx.size(); ++idxq) {
        double tmpval = modesum_sv(mode_coefficients, mode_all, inp_model,
                                   idxlayer, idxpl, idxq);
        double xscale = delta * qx[idxq] + xpl;
        tmpval -= modecomp.v(idxlayer)(xscale);
        tmpsum += delta * tmpval * tmpval * qw[idxq];
      }
      totint += tmpsum;
    };
  };
  return std::sqrt(totint);
};

auto
traction_singlemode(ModeCoupling::mineos_eigenfunction_continuous &mode1,
                    EarthModels::ModelInput<double, int> &inp_model) {
  int idxlouter = inp_model.NumberOfLayers() - 1;
  double r = inp_model.OuterRadius();
  double kappa =
      inp_model.Density(idxlouter)(r) *
      (inp_model.VP(idxlouter)(r) * inp_model.VP(idxlouter)(r) -
       4.0 / 3.0 * inp_model.VS(idxlouter)(r) * inp_model.VS(idxlouter)(r));
  double mu = inp_model.Density(idxlouter)(r) * inp_model.VS(idxlouter)(r) *
              inp_model.VS(idxlouter)(r);

  double tractionval = 0.0;

  tractionval = (kappa + 4.0 / 3.0 * mu) * mode1.up(idxlouter)(r) +
                (kappa - 2.0 / 3.0 * mu) *
                    (2.0 * mode1.u(idxlouter)(r) -
                     std::sqrt(2.0) * mode1.v(idxlouter)(r)) /
                    r * inp_model.OuterRadius();
  return tractionval;
};

auto
traction_singlemode2(ModeCoupling::mineos_eigenfunction_continuous &mode1,
                     EarthModels::ModelInput<double, int> &inp_model, int idx,
                     double r) {
  // int idxlouter = inp_model.NumberOfLayers() - 1;
  // double r = inp_model.OuterRadius();
  double kappa = inp_model.Density(idx)(r) *
                 (inp_model.VP(idx)(r) * inp_model.VP(idx)(r) -
                  4.0 / 3.0 * inp_model.VS(idx)(r) * inp_model.VS(idx)(r));
  double mu =
      inp_model.Density(idx)(r) * inp_model.VS(idx)(r) * inp_model.VS(idx)(r);

  double tractionval = 0.0;

  tractionval = (kappa + 4.0 / 3.0 * mu) * mode1.up(idx)(r) +
                (kappa - 2.0 / 3.0 * mu) *
                    (2.0 * mode1.u(idx)(r) - std::sqrt(2.0) * mode1.v(idx)(r)) /
                    r * inp_model.OuterRadius();
  return tractionval;
};

auto
traction_basis_sum(Eigen::VectorXd &mode_coefficients,
                   ModeCoupling::modesummary &mode_all,
                   EarthModels::ModelInput<double, int> &inp_model) {

  int idxlouter = inp_model.NumberOfLayers() - 1;
  int idxplouter = inp_model.LayerRadii(idxlouter).size() - 2;
  int idxqouter = mode_all.npoints() - 1;
  double r = inp_model.OuterRadius();
  double kappa =
      inp_model.Density(idxlouter)(r) *
      (inp_model.VP(idxlouter)(r) * inp_model.VP(idxlouter)(r) -
       4.0 / 3.0 * inp_model.VS(idxlouter)(r) * inp_model.VS(idxlouter)(r));
  double mu = inp_model.Density(idxlouter)(r) * inp_model.VS(idxlouter)(r) *
              inp_model.VS(idxlouter)(r);

  double tractionval = 0.0;
  tractionval += (kappa + 4.0 / 3.0 * mu) *
                 modesum_sup(mode_coefficients, mode_all, inp_model, idxlouter,
                             idxplouter, idxqouter);
  tractionval +=
      (kappa - 2.0 / 3.0 * mu) *
      (2.0 * modesum_s(mode_coefficients, mode_all, inp_model, idxlouter,
                       idxplouter, idxqouter) -
       std::sqrt(2.0) * modesum_sv(mode_coefficients, mode_all, inp_model,
                                   idxlouter, idxplouter, idxqouter)) /
      r * inp_model.OuterRadius();
  return tractionval;
};

auto
h10norm(ModeCoupling::mineos_eigenfunction_continuous &mode1,
        EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = q.Points();
  std::vector<double> qw = q.Weights();

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
      auto fun = [&inp_model, &mode1, &idxlayer, &vec_radii, &idxpl](double x) {
        double x1 = vec_radii[idxpl];
        double x2 = vec_radii[idxpl + 1];

        double delta = (x2 - x1) / 2.0;
        double xpl = (x1 + x2) / 2.0;
        double xscale = delta * x + xpl;
        double tmp1 = mode1.u(idxlayer)(xscale);
        double tmp2 = mode1.up(idxlayer)(xscale);
        return delta * (tmp1 * tmp1 + tmp2 * tmp2);
      };
      totint += q.Integrate(fun);
    };
  };
  return std::sqrt(totint);
};

auto
h10normdifference_s(Eigen::VectorXd &mode_coefficients,
                    ModeCoupling::modesummary &mode_all,
                    ModeCoupling::mineos_eigenfunction_continuous &modecomp,
                    EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = mode_all.q().Points();
  std::vector<double> qw = mode_all.q().Weights();
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
      for (int idxq = 0; idxq < qx.size(); ++idxq) {
        double xscale = delta * qx[idxq] + xpl;
        double tmpval1 = modesum_s(mode_coefficients, mode_all, inp_model,
                                   idxlayer, idxpl, idxq);
        tmpval1 -= modecomp.u(idxlayer)(xscale);
        double tmpval2 = modesum_sup(mode_coefficients, mode_all, inp_model,
                                     idxlayer, idxpl, idxq);
        tmpval2 -= modecomp.up(idxlayer)(xscale);
        tmpsum += delta * qw[idxq] * (tmpval1 * tmpval1 + tmpval2 * tmpval2);
      }
      totint += tmpsum;
    };
  };
  return std::sqrt(totint);
};

auto
h10normdifference_s2(Eigen::VectorXd &mode_coefficients,
                     ModeCoupling::modesummary2 &modes_all,
                     ModeCoupling::mineos_eigenfunction_continuous &modecomp,
                     EarthModels::ModelInput<double, int> &inp_model) {
  // Try to use integrate in gauss quad
  // auto q = GaussQuad::GaussLobattoLegendreQuadrature1D<double>(6);
  std::vector<double> qx = modes_all.q().Points();
  std::vector<double> qw = modes_all.q().Weights();
  double totint = 0.0;
  double densitynorm = 5515.0;
  double radiusnorm = inp_model.OuterRadius() * inp_model.LengthNorm();
  //    std::cout << radiusnorm << "\n\n";
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double velocitynorm = radiusnorm / std::sqrt(pi_db * bigg_db * densitynorm);
  double frequencynorm = velocitynorm / radiusnorm;

  for (int idxe = 0; idxe < modes_all.fullmesh().NE(); ++idxe) {
    int idxlayer = modes_all.fullmesh().LayerNumber(idxe);
    double x1 = modes_all.fullmesh().ELR(idxe);
    double x2 = modes_all.fullmesh().EUR(idxe);
    // if ((idx1 + idx2) == 0) {
    //    std::cout << idxlayer << " " << x1 << " " << x2 << "\n";
    // }
    double delta = (x2 - x1) / 2.0;
    double xpl = (x1 + x2) / 2.0;
    double tmpsum = 0.0;
    for (int idxq = 0; idxq < qx.size(); ++idxq) {

      double xscale = delta * qx[idxq] + xpl;
      double tmpval =
          modesum_s2(mode_coefficients, modes_all, inp_model, idxe, idxq);
      tmpval -= modecomp.u(idxlayer)(xscale);

      double tmpval2 =
          modesum_sv2(mode_coefficients, modes_all, inp_model, idxe, idxq);
      tmpval2 -= modecomp.v(idxlayer)(xscale);

      double tmpval3 =
          modesum_sup2(mode_coefficients, modes_all, inp_model, idxe, idxq);
      tmpval3 -= modecomp.up(idxlayer)(xscale);

      double tmpval4 =
          modesum_svp2(mode_coefficients, modes_all, inp_model, idxe, idxq);
      tmpval4 -= modecomp.vp(idxlayer)(xscale);

      tmpsum += delta * qw[idxq] *
                (tmpval * tmpval + tmpval2 * tmpval2 + tmpval3 * tmpval3 +
                 tmpval4 * tmpval4);
    }
    totint += tmpsum;
  };
  return std::sqrt(totint);
};
