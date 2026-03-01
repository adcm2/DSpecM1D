#ifndef TOROIDAL_GAL_CLASS_GUARD_H
#define TOROIDAL_GAL_CLASS_GUARD_H
#include <iostream>
#include <cmath>
#include <functional>
#include <GaussQuad/All>
#include <Interpolation/All>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymGEigsSolver.h>
#include <PlanetaryModel/All>
#include <Spectra/SymGEigsShiftSolver.h>

using namespace Spectra;

namespace Toroidal {

class gal_gen {
private:
  Eigen::VectorXcd _evalues;
  Eigen::MatrixXcd _evectors, mat_inertia, mat_ke;
  std::vector<std::vector<double>> _nodes;
  GaussQuad::Quadrature1D<double> _q;
  std::size_t _numbasis, _numeig;
  std::size_t index_evec(std::size_t i, std::size_t k) {
    return (_q.N() - 1) * i + k;
  };
  Eigen::VectorXcd vec_eig;
  Eigen::MatrixXcd mat_eig;

  using VCD = std::vector<std::complex<double>>;
  using VVCD = std::vector<VCD>;
  using VVVCD = std::vector<VVCD>;
  using VD = std::vector<double>;
  using VVD = std::vector<VD>;
  using VVVD = std::vector<VVD>;

  VVVD val_ef, val_efd;
  VVVCD val_efc, val_efdc;
  //   VVVD val_full,val_;
  double densitynorm = 5515.0;
  double pi_db = 3.14159265358979;
  double bigg_db = 6.6723 * std::pow(10.0, -11.0);
  double frequencynorm = std::sqrt(pi_db * bigg_db * densitynorm);
  double normint;
  double _freq_norm;

public:
  gal_gen() {};
  template <class model1d>
  gal_gen(spectral_element_planet &, const model1d &, const model1d &, bool);
  //   gal_gen(const spectral_element &, const model_string &, const
  //   model_string &,
  //           bool);

  // output
  auto Modes() const { return val_efc; };
  auto Mode_Derivatives() const { return val_efdc; };
  VVVCD &Modes_Ref() { return val_efc; };

  // traction
  template <class model1d>
  VVVCD Traction(spectral_element_planet &inp_basis,
                 const model1d &inp_model) const {
    auto _mesh = inp_basis.mesh();
    std::size_t nume = _mesh.NE();
    std::size_t numq = _mesh.NN();
    std::cout << "\nne: " << nume << ", numq: " << numq
              << ", numb: " << _numbasis << "\n\n";
    VVVCD vec_traction = VVVCD(nume, VVCD(numq, VCD(_numbasis, 0.0)));

    for (int idxe = 0; idxe < nume; ++idxe) {
      int qmin = (idxe == 0);
      for (int idxq = qmin; idxq < numq; ++idxq) {
        double cr = _mesh.NodeRadius(idxe, idxq);
        for (int i = 0; i < _numbasis; ++i) {
          vec_traction[idxe][idxq][i] =
              inp_model.L(_mesh.LayerNumber(idxe))(cr) *
              (val_efdc[idxe][idxq][i] - val_efc[idxe][idxq][i] / cr);
        }
      }
    }
    return vec_traction;
  };

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  auto specval(double fv, double ep, Eigen::VectorXcd &vec_f,
               Eigen::VectorXcd &vec_val) {
    std::complex<double> wv = (fv, -ep);
    Eigen::MatrixXcd mat_w = mat_ke - wv * wv * mat_inertia;
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(mat_w);
    Eigen::VectorXcd vec_x = lu.solve(vec_f);
    auto tmp = vec_val.adjoint() * vec_x;
    return tmp;
  }

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  Eigen::VectorXd efrequencies() {
    Eigen::VectorXd ret_eig(vec_eig.rows());
    for (int i = 0; i < vec_eig.rows(); ++i) {
      ret_eig(i) = std::sqrt(std::abs(vec_eig(i)));
    }
    return ret_eig;
  };
  Eigen::VectorXcd evalues() { return vec_eig; };
  Eigen::MatrixXcd evectors() { return mat_eig; };
};

template <class model1d>
gal_gen::gal_gen(spectral_element_planet &inp_basis, const model1d &mod_init,
                 const model1d &mod_new, bool to_aug)
    : _evalues{inp_basis.evalues()}, _evectors{inp_basis.evectors_all()},
      _freq_norm{1.0 / mod_new.TimeNorm()} {

  // std::cout << "Frequency: " << _freq_norm << ", 2: " << frequencynorm
  //           << ", density: " << mod_new.DensityNorm() << ", 2: " <<
  //           densitynorm
  //           << " \n";
  normint =
      _freq_norm / frequencynorm * sqrt(mod_new.DensityNorm() / densitynorm);
  // std::cout << ", normint: " << normint << "\n";
  // check that the new model has the same discontinuities
  if (mod_init.NumberOfLayers() == mod_new.NumberOfLayers()) {
    for (int idx = 0; idx < mod_init.NumberOfLayers(); ++idx) {
      if ((mod_init.UpperRadius(idx) - mod_new.UpperRadius(idx)) >
          std::pow(10.0, -12.0)) {
        std::cout << "Models have different discontinuities!\n";
      }
    }
  } else {
    std::cout << "Models have different number of layers!\n";
  }
  assert(inp_basis.EigenCalc() && "Eigenmodes not calculated.");

  auto q = inp_basis.q();
  // checked that the models have the same discontinuities
  // find size of matrix
  int el = inp_basis.el();
  int eu = inp_basis.eu();
  std::size_t matlen = inp_basis.idx_submesh(eu - 1, q.N());
  // std::cout << "\nSize: " << matlen << "\n";

  // size is different if we augment or not
  // int _numbasis;
  if (to_aug) {
    _numbasis = inp_basis.NumberOfModes() + inp_basis.NumberOfAugment();
    std::cout << "\nSize of basis 1: " << _numbasis << ", aug: " << to_aug
              << "\n\n";
  } else {
    _numbasis = inp_basis.NumberOfModes();
    std::cout << "\nSize of basis 2: " << _numbasis << ", aug: " << to_aug
              << "\n\n";
  }

  // resize the matrices
  mat_ke.resize(_numbasis, _numbasis);
  mat_inertia.resize(_numbasis, _numbasis);

  // modes and its derivative
  auto evectors = inp_basis.evectors_std_ref();
  auto evectors_deriv = inp_basis.evectors_deriv_ref();
  auto aug_val = inp_basis.aug_std_ref();
  auto aug_deriv = inp_basis.aug_deriv_ref();

  // mesh
  auto _mesh = inp_basis.mesh();

  // fill out val_ef and val_efd:
  val_ef = VVVD(_mesh.NE(), VVD(q.N(), VD(_numbasis, 0.0)));
  val_efd = val_ef;
  val_efc = VVVCD(_mesh.NE(), VVCD(q.N(), VCD(_numbasis, 0.0)));
  val_efdc = val_efc;

  // std::cout << "Check 1\n";
  if (to_aug) {
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      for (int idxq = 0; idxq < q.N(); ++idxq) {
        for (int i = 0; i < inp_basis.NumberOfModes(); ++i) {
          val_ef[idxe][idxq][i] = evectors[idxe][idxq][i];
          val_efd[idxe][idxq][i] = evectors_deriv[idxe][idxq][i];
        }
        for (int i = inp_basis.NumberOfModes(); i < _numbasis; ++i) {
          val_ef[idxe][idxq][i] =
              aug_val[idxe][idxq][i - inp_basis.NumberOfModes()];
          val_efd[idxe][idxq][i] =
              aug_deriv[idxe][idxq][i - inp_basis.NumberOfModes()];
        }
      }
    }
  } else {
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      for (int idxq = 0; idxq < q.N(); ++idxq) {
        for (int i = 0; i < inp_basis.NumberOfModes(); ++i) {
          val_ef[idxe][idxq][i] = evectors[idxe][idxq][i];
          val_efd[idxe][idxq][i] = evectors_deriv[idxe][idxq][i];
        }
      }
    }
  }

  // inertia matrix
  for (int i = 0; i < _numbasis; ++i) {
    for (int j = 0; j < i + 1; ++j) {
      double tmp = 0.0;
      // integrate over each element
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        double elem_width = _mesh.EW(idxe);
        double d_val = 2.0 / elem_width;
        int laynum = _mesh.LayerNumber(idxe);
        double tmp1 = 0.0;
        for (int idxq = 0; idxq < q.N(); ++idxq) {
          double cr = _mesh.NodeRadius(idxe, idxq);
          double tmp2 = val_ef[idxe][idxq][i] * val_ef[idxe][idxq][j];
          tmp2 *= mod_new.Density(laynum)(cr) * cr * cr;
          tmp1 += q.W(idxq) * tmp2;
        }
        tmp += tmp1 * elem_width * 0.5;
      }
      if (i != j) {
        mat_inertia(i, j) = tmp;
        mat_inertia(j, i) = tmp;
      } else {
        mat_inertia(i, j) = tmp;
      }
    }
  }

  // ke matrix
  int _k2 = inp_basis.k2();
  for (int i = 0; i < _numbasis; ++i) {
    for (int j = 0; j < i + 1; ++j) {
      double tmp = 0.0;
      // integrate over each element
      for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
        double elem_width = _mesh.EW(idxe);
        int laynum = _mesh.LayerNumber(idxe);
        double tmp1 = 0.0;
        for (int idxq = 0; idxq < q.N(); ++idxq) {
          double cr = _mesh.NodeRadius(idxe, idxq);

          // first integral
          double tmp2 = cr * val_efd[idxe][idxq][i] - val_ef[idxe][idxq][i];
          double tmp3 = cr * val_efd[idxe][idxq][j] - val_ef[idxe][idxq][j];
          double tmp4 = tmp2 * tmp3 * mod_new.L(laynum)(cr);

          // second integral
          double tmp5 =
              (_k2 - 2.0) * val_ef[idxe][idxq][i] * val_ef[idxe][idxq][j];
          tmp5 *= mod_new.N(laynum)(cr);

          // sum
          tmp1 += q.W(idxq) * (tmp4 + tmp5);
        }
        tmp += tmp1 * elem_width * 0.5;
      }
      if (i != j) {
        mat_ke(i, j) = tmp;
        mat_ke(j, i) = tmp;
      } else {
        mat_ke(i, j) = tmp;
      }
    }
  }

  // std::cout << "Check 2\n";

  // setting up the problem, we need to find the inverse of the inertia matrix
  // and pre-multiply it through against the kinetic energy matrix
  Eigen::FullPivLU<Eigen::MatrixXcd> lu(mat_inertia);
  Eigen::MatrixXcd mat_decomp = lu.solve(mat_ke);

  // eigendecompose
  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(mat_decomp);
  vec_eig = es.eigenvalues();
  mat_eig = es.eigenvectors();

  // normalisation
  Eigen::MatrixXcd mat_norm = mat_eig.transpose() * mat_inertia * mat_eig;
  for (int i = 0; i < mat_norm.rows(); ++i) {
    std::cout << std::setprecision(2) << i << " " << mat_norm(i, i) << "\n";
  }
  std::cout << "\n" << "normint: " << normint << "\n\n";
  // normalise
  for (int i = 0; i < mat_eig.rows(); ++i) {
    for (int j = 0; j < mat_eig.cols(); ++j) {
      mat_eig(i, j) *= 1.0 / std::sqrt(mat_norm(j, j).real());
    }
  }

  // normalise:
  for (int i = 0; i < mat_eig.rows(); ++i) {
    for (int j = 0; j < mat_eig.cols(); ++j) {
      mat_eig(i, j) *= 1.0 / (normint * sqrt(vec_eig(j).real()));
    }
  }

  // calculate eigenvectors

  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < q.N(); ++idxq) {
      for (int i = 0; i < _numbasis; ++i) {
        std::complex<double> tmp = 0.0, tmp1 = 0.0;
        for (int j = 0; j < _numbasis; ++j) {
          tmp += mat_eig(j, i) * val_ef[idxe][idxq][j];
          tmp1 += mat_eig(j, i) * val_efd[idxe][idxq][j];
        }
        val_efc[idxe][idxq][i] = tmp;
        val_efdc[idxe][idxq][i] = tmp1;
      }
    }
  }
};

}   // namespace Toroidal

#endif