#ifndef MESH_MODEL_H
#define MESH_MODEL_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <EarthMesh/All>

class MeshModel {
public:
  template <class model1d> MeshModel(EarthMesh::RadialMesh &, model1d &);
  MeshModel() {};
  //   ~MeshModel();

  auto Density(int idxe, int idxn) { return _vec_density[idxe][idxn]; };
  auto VPV(int idxe, int idxn) { return _vec_vpv[idxe][idxn]; };
  auto VSV(int idxe, int idxn) { return _vec_vsv[idxe][idxn]; };
  auto VPH(int idxe, int idxn) { return _vec_vph[idxe][idxn]; };
  auto VSH(int idxe, int idxn) { return _vec_vsh[idxe][idxn]; };
  auto VP(int idxe, int idxn) { return _vec_vp[idxe][idxn]; };
  auto VS(int idxe, int idxn) { return _vec_vs[idxe][idxn]; };
  auto QKappa(int idxe, int idxn) { return _vec_qkappa[idxe][idxn]; };
  auto QMu(int idxe, int idxn) { return _vec_qmu[idxe][idxn]; };
  auto Eta(int idxe, int idxn) { return _vec_eta[idxe][idxn]; };
  auto A(int idxe, int idxn) { return _vec_A[idxe][idxn]; };
  auto N(int idxe, int idxn) { return _vec_N[idxe][idxn]; };
  auto L(int idxe, int idxn) { return _vec_L[idxe][idxn]; };
  auto C(int idxe, int idxn) { return _vec_C[idxe][idxn]; };
  auto F(int idxe, int idxn) { return _vec_F[idxe][idxn]; };
  auto Kappa(int idxe, int idxn) { return _vec_Kappa[idxe][idxn]; };
  auto Mu(int idxe, int idxn) { return _vec_Mu[idxe][idxn]; };
  auto PSlow(int idxe, int idxn) { return _vec_pslow[idxe][idxn]; };
  auto SSlow(int idxe, int idxn) { return _vec_sslow[idxe][idxn]; };

  // attenuation equivalent parameters
  auto A_atten(int idxe, int idxn) { return _vec_A_atten[idxe][idxn]; };
  auto N_atten(int idxe, int idxn) { return _vec_N_atten[idxe][idxn]; };
  auto L_atten(int idxe, int idxn) { return _vec_L_atten[idxe][idxn]; };
  auto C_atten(int idxe, int idxn) { return _vec_C_atten[idxe][idxn]; };
  auto F_atten(int idxe, int idxn) { return _vec_F_atten[idxe][idxn]; };

  // gravity
  auto Gravity(int idxe, int idxn) { return _vec_gravity[idxe][idxn]; };

private:
  using VECTOR = std::vector<double>;
  using VVECTOR = std::vector<std::vector<double>>;
  double fov3 = 4.0 / 3.0;
  double tov3 = 2.0 / 3.0;
  VVECTOR _vec_density;
  VVECTOR _vec_vpv;
  VVECTOR _vec_vsv;
  VVECTOR _vec_vph;
  VVECTOR _vec_vsh;
  VVECTOR _vec_vp;
  VVECTOR _vec_vs;
  VVECTOR _vec_qkappa;
  VVECTOR _vec_qmu;
  VVECTOR _vec_eta;
  VVECTOR _vec_A;
  VVECTOR _vec_N;
  VVECTOR _vec_L;
  VVECTOR _vec_C;
  VVECTOR _vec_F;
  VVECTOR _vec_Kappa;
  VVECTOR _vec_Mu;
  VVECTOR _vec_A_atten;
  VVECTOR _vec_N_atten;
  VVECTOR _vec_L_atten;
  VVECTOR _vec_C_atten;
  VVECTOR _vec_F_atten;
  VVECTOR _vec_gravity;
  VVECTOR _vec_pslow;
  VVECTOR _vec_sslow;
};

template <class model1d>
MeshModel::MeshModel(EarthMesh::RadialMesh &mesh, model1d &inp_model) {
  int NE = mesh.NE();
  int NQ = mesh.NN();

  _vec_density = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_vpv = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_vsv = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_vph = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_vsh = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_vp = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_vs = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_qkappa = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_qmu = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_eta = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_A = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_N = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_L = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_C = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_F = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_Kappa = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_Mu = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_A_atten = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_N_atten = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_L_atten = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_C_atten = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_F_atten = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_pslow = VVECTOR(NE, VECTOR(NQ, 0.0));
  _vec_sslow = VVECTOR(NE, VECTOR(NQ, 0.0));

  for (int idxe = 0; idxe < NE; ++idxe) {
    auto laynum = mesh.LayerNumber(idxe);
    for (int idxn = 0; idxn < NQ; ++idxn) {
      double r = mesh.NodeRadius(idxe, idxn);
      _vec_density[idxe][idxn] = inp_model.Density(laynum)(r);
      _vec_vpv[idxe][idxn] = inp_model.VPV(laynum)(r);
      _vec_vsv[idxe][idxn] = inp_model.VSV(laynum)(r);
      _vec_vph[idxe][idxn] = inp_model.VPH(laynum)(r);
      _vec_vsh[idxe][idxn] = inp_model.VSH(laynum)(r);
      _vec_vp[idxe][idxn] = inp_model.VP(laynum)(r);
      _vec_vs[idxe][idxn] = inp_model.VS(laynum)(r);
      _vec_qkappa[idxe][idxn] = inp_model.QKappa(laynum)(r);
      _vec_qmu[idxe][idxn] = inp_model.QMu(laynum)(r);
      _vec_eta[idxe][idxn] = inp_model.Eta(laynum)(r);
      _vec_A[idxe][idxn] = inp_model.A(laynum)(r);
      _vec_N[idxe][idxn] = inp_model.N(laynum)(r);
      _vec_L[idxe][idxn] = inp_model.L(laynum)(r);
      _vec_C[idxe][idxn] = inp_model.C(laynum)(r);
      _vec_F[idxe][idxn] = inp_model.F(laynum)(r);
      _vec_Kappa[idxe][idxn] = inp_model.Kappa(laynum)(r);
      _vec_Mu[idxe][idxn] = inp_model.Mu(laynum)(r);
      _vec_pslow[idxe][idxn] =
          1.0 / (_vec_vp[idxe][idxn] * _vec_vp[idxe][idxn]);
      _vec_sslow[idxe][idxn] =
          1.0 / (_vec_vs[idxe][idxn] * _vec_vs[idxe][idxn]);

      // equivalent attenuation parameters
      auto ratkappa = _vec_Kappa[idxe][idxn] / _vec_qkappa[idxe][idxn];
      auto ratmu = _vec_Mu[idxe][idxn] / _vec_qmu[idxe][idxn];
      if (_vec_qmu[idxe][idxn] != 0.0) {
        _vec_A_atten[idxe][idxn] += fov3 * ratmu;
        _vec_C_atten[idxe][idxn] += fov3 * ratmu;
        _vec_N_atten[idxe][idxn] += ratmu;
        _vec_L_atten[idxe][idxn] += ratmu;
        _vec_F_atten[idxe][idxn] += -tov3 * ratmu;
      }
      if (_vec_qkappa[idxe][idxn] != 0.0) {
        _vec_A_atten[idxe][idxn] += ratkappa;
        _vec_C_atten[idxe][idxn] += ratkappa;
        _vec_F_atten[idxe][idxn] += ratkappa;
      }
    }
  }

  _vec_gravity = VVECTOR(NE, VECTOR(NQ, 0.0));
  auto bigg_db =
      6.67230 * std::pow(10.0, -11.0) / inp_model.GravitationalConstant();
  auto pi_db = 3.14159265358979323846;
  auto q = mesh.GLL();
  // compute gravity at all nodes
  for (int idxe = 0; idxe < NE; ++idxe) {
    int laynum = mesh.LayerNumber(idxe);
    if (idxe != 0) {
      _vec_gravity[idxe][0] = _vec_gravity[idxe - 1][NQ - 1];
    }
    int idxlow = 1;
    for (int idxn = idxlow; idxn < NQ; ++idxn) {
      double urad = mesh.NodeRadius(idxe, idxn);
      double lrad = mesh.NodeRadius(idxe, idxn - 1);
      double tmp = 0.0;
      double ewidth2 = 0.5 * (urad - lrad);
      double crad2 = 0.5 * (urad + lrad);
      for (int idxi = 0; idxi < NQ; ++idxi) {
        double cradi = ewidth2 * q.X(idxi) + crad2;
        tmp += q.W(idxi) * _vec_density[idxe][idxi] * cradi * cradi;
      }
      tmp *= ewidth2;
      _vec_gravity[idxe][idxn] = _vec_gravity[idxe][idxn - 1] + tmp;
    }
  }

  // divide by r^2:
  for (int idxe = 0; idxe < NE; ++idxe) {
    int idxlow = (idxe == 0);
    for (int idxn = idxlow; idxn < NQ; ++idxn) {
      double crad = mesh.NodeRadius(idxe, idxn);
      _vec_gravity[idxe][idxn] *= 4.0 * pi_db * bigg_db / (crad * crad);
    }
  }
}

#endif   // MESH_MODEL_H